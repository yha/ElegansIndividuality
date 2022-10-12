using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

##

using CSV, DataFrames
using OrderedCollections
using NamedArrays, BlockArrays
using StatsBase
using IterTools
using Random
using LinearAlgebra
using MultivariateStats
using RLEVectors
using JLD2

include("src/funcs.jl")

const STAGE_NAMES = ["L1", "L2", "L3", "L4", "A"]
const N_STAGES = length(STAGE_NAMES)


## Load data

nbins = 10

well_conds_df = CSV.read("data/well-list.csv", DataFrame)
roam_frac_df = CSV.read("data/roam-frac-$nbins.csv", DataFrame)

wells = [(; row.experiment, row.well) for row in eachrow(well_conds_df)]
@assert wells == [(; row.experiment, row.well) for row in eachrow(roam_frac_df)]

roam_frac_mat = Matrix(roam_frac_df[:,3:end]) :: Matrix{Float64}

@assert size(roam_frac_mat, 2) == N_STAGES * nbins

## Prepare roaming fraction matrices per condition

conditions = sort!(unique(zip(well_conds_df.strain, well_conds_df.DS)))
cond2wells = OrderedDict(
    (strain, ds) => [(; r.experiment, r.well) 
                            for r in eachrow(well_conds_df) if (r.strain, r.DS) == (strain, ds)]
    for (strain, ds) in conditions)


cond2well_i = OrderedDict( c => Vector{Int}(indexin(cond2wells[c], wells)) for c in conditions )

cond2roam_frac_mat = OrderedDict( c => view(roam_frac_mat,i,:) for (c,i) in cond2well_i )

function condition_roam_mat( condition, stage; nbins, cond2roam_frac_mat )
	m = cond2roam_frac_mat[condition]
	@assert size(m,2) == 5nbins
	m[:, nbins*(stage-1)+1:nbins*stage]
end

# Retrieve roaming fraction matrix for given conditions and stages, 
# as a block matrix with named blocks per condition and stage
function get_roam_mat_and_wells( conditions, stages; 
							  nbins = nbins, cond2roam_frac_mat, cond2wells )
	m = [condition_roam_mat(c, stage; nbins, cond2roam_frac_mat) 
		 for c in conditions, stage in stages]
   	condnames = join.(conditions, " ")
	stagenames = STAGE_NAMES[stages]
	
    fracs   = NamedArray( m,   (condnames, stagenames), ("Conditions", "Stages") )
	wells 	= NamedArray( [cond2wells[c] for c in conditions], condnames, "Conditions" )
	mortar(fracs), mortar(wells)
end


## Compute and save output

# Compute everything for given conditions and stages
function pipeline(conditions, stages; 
					cond2roam_frac_mat, cond2wells, nshuffles, n_pair_row_shuffles)
	cond_str = join(("$strain:$(ds)DS" for (strain, ds) in conditions), ", ")
	@info "pipeline for $cond_str..." stages
	@info "...roaming fraction and ranks"
	roamfrac, wells = get_roam_mat_and_wells( conditions, stages; 
								cond2roam_frac_mat, cond2wells )
	ranks = rank_each_experiment( roamfrac, wells )

	@info "...bin shuffles"
	@time shuffled_roamfrac = [shuffle_per_experiment(roamfrac, wells) for _ = 1:nshuffles]
	@time shuffled_ranks = [rank_each_experiment(m, wells) for m in shuffled_roamfrac]

	@info "...consistency"
	consistency = consistency_index.(eachrow(ranks))
	shuffled_consistency = [consistency_index.(eachrow(sh)) 
	for sh in shuffled_ranks]

	@info "...condition pair individual shuffles"
	@time pair_shuffles = make_cond_pair_shuffles(ranks; 
				nshuffles = n_pair_row_shuffles)

	@info "...PCA"
	pca = multi_cond_pca(ranks)
	coeffs, P, μ = positive_transform(pca, ranks')

	(; wells, roamfrac, ranks, shuffled_roamfrac, shuffled_ranks, consistency, shuffled_consistency, pair_shuffles, pca, coeffs, P, μ)
end

function make_rank_and_coeff_csvs(prefix, conditions, wells, ranks, coeffs)
	row_conds = RLEVector(conditions, blocklasts(axes(ranks,1)))
	common_columns = hcat(DataFrame(wells), 
					 	  DataFrame([(;strain=r[1], DS=r[2]) for r in row_conds]))

	ranks_df =  hcat(common_columns, 
				    DataFrame(Tables.table(ranks; 
				   			  header = "bin " .* string.(1:size(ranks,2)))))

	scores_df = hcat(common_columns, 
					 DataFrame(coeffs', "PC " .* string.(1:size(coeffs,1))))

	@info "Saving ranks and score CSV files for $prefix"
	CSV.write("$prefix ranks.csv",  ranks_df)
	CSV.write("$prefix scores.csv", scores_df)
end

function make_and_save_data(conditions, prefix; kwargs... )
	r = pipeline( conditions, 1:5; kwargs... );
	make_rank_and_coeff_csvs(prefix, conditions, r.wells, r.ranks, r.coeffs)
	fname = "$prefix data.jld2"
	@info "Saving data to \"$fname\""
	jldsave(fname; r...)
end

nshuffles, n_pair_row_shuffles = 100, 100
outdir = "output"
mkpath(outdir)

let N2_conds = [(strain,ds) for (strain,ds) in conditions if strain == "N2"]
	make_and_save_data(N2_conds, "$outdir/N2"; 
		cond2roam_frac_mat, cond2wells, nshuffles, n_pair_row_shuffles)
end

let conds = [(strain,ds) for (strain,ds) in conditions if strain ∈ ("N2", "tph-1", "cat-2")]
	make_and_save_data(conds, "$outdir/N2, tph-1, cat-2"; 
		cond2roam_frac_mat, cond2wells, nshuffles, n_pair_row_shuffles)
end

