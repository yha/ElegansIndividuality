## Utils for named-block arrays with conditions in block-rows and stages in block-columns

with_named_blocks(block_array) = mortar(NamedArray(identity.(blocks(block_array))))
name_blocks(block_array, names) = mortar(NamedArray(identity.(blocks(block_array)), names))
function name_blocks(block_array, names, d)
	a = with_named_blocks(block_array)
	setnames!(blocks(a), names, d)
	a
end

function each_experiment(mat, wells)
	@assert size(mat,1) == length(wells)
	_, lens = experiment_rle(wells)
	ex_ends = cumsum(lens)
	(mat[i+1:j,:] for (i,j) in IterTools.partition( [0; ex_ends], 2, 1 ))
end

function experiment_rle(wells)
	well_exs = [w.experiment for w in vec(wells)]
	exs, lens = rle(well_exs)
	# each experiments should appear in one continuous block
	@assert unique(well_exs) == exs
	exs, lens
end

mimic_blocks(a, ba) = name_blocks(PseudoBlockArray(a, axes(ba)), names(blocks(ba)))

shuffle_each_bin(means) = mapslices(shuffle, means; dims=1)
shuffle_each_bin_and_fake_blocks(means) = mimic_blocks(shuffle_each_bin(means), means)

shuffle_per_condition(means) = mortar(NamedArray(reduce(vcat, blocks(shuffle_each_bin(c)) for c in eachblockrow(means)), names(blocks(means))))

function shuffle_per_experiment(means, wells)
    mat = reduce(vcat, shuffle_each_bin(c) for c in each_experiment(means, wells))
    mimic_blocks(mat, means)
end

shuffled(blockmat, i) = mimic_blocks(blockmat[i], blockmat)

## Ranking

rank_rescale(r,n) = 2/n * (r-1/2) - 1
rescale_rank_vec(v) = rank_rescale.(v, length(v))
rank_each_col(x) = reduce(hcat, rescale_rank_vec(tiedrank(c)) for c in eachcol(x))

function rank_each_experiment( mat, wells )
	# Ranking applied to each experiment×stage combination.
	#    This is done in each block (condition×stage) separately 
	#    to preserve the block structure.
	normed_mblocks = [reduce(vcat, 
			map(rank_each_col, each_experiment(mat[i,j], wells[i]))) 	
		for i in blockaxes(mat,1), j in blockaxes(mat,2)]	
	
	mortar(NamedArray(normed_mblocks, names(blocks(mat))))
end

## Shuffling

blockranges(a,dim) = range.(blockfirsts(axes(a,dim)), blocklasts(axes(a,dim)))
reshuffle(v,w) = let k = shuffle(vcat(v,w))
    k[1:length(v)], k[length(v)+1:end]
end
function make_row_shuffles(mat, block_i, block_j, n)
    r1, r2 = blockranges(mat, 1)[[block_i, block_j]]
    map(1:n) do _
        rs1, rs2 = reshuffle(r1, r2)
        #mat[rs1,:], mat[rs2,:]
        view(mat,rs1,:), view(mat,rs2,:)
    end
end

function make_cond_pair_shuffles(blockmat; nshuffles)
	nconds = length(blockaxes(blockmat, 1))
	# compute shuffles for pairs j <= i
	# `sh_ik, sh_jk = S[i][j][k]` is the k-th shuffle,
	# with `sh_ik` standing for condition `i`, and `sh_jk` standing for condition `j`
	S = [[make_row_shuffles(blockmat, i, j, nshuffles) for j=1:i] 
		for i=1:nconds]

	# For i<j, the (i,j) shuffles are the (j,i) with each pair reversed
	return get_pair_shuffles(i,j) = if j <= i
		S[i][j]
	else
		[sh[[2,1]] for sh in S[j][i]]
	end
end

## Weighted multi-condition PCA

function weighted_pca( X, weights; kwargs... )
    X isa Adjoint && (X = Matrix(X)) # workaround for Statistics issue #82
    μ, Σ = mean_and_cov(X, FrequencyWeights(weights))
    pcacov( Σ, μ[1,:]; kwargs... )
end

function multi_cond_pca(M, cond_ends, cond_weights; pca_kwargs...)
    cond_lens = diff([0; cond_ends])
    row_weights = inverse_rle(cond_weights ./ cond_lens, cond_lens)
    weighted_pca(M, row_weights; pca_kwargs...)
end

function multi_cond_pca(M::AbstractBlockMatrix, 
                        cond_weights = ones(length(blockaxes(M,1))); 
                        pca_kwargs...)
    multi_cond_pca(Matrix(M), blocklasts(axes(M,1)), cond_weights; pca_kwargs...)
end

function positive_projection(pca; f=sum)
    proj = projection(pca)
    signs = sign.(f.(eachcol(proj)))
    proj .* signs'
end

function positive_transform(pca, m; f=sum)
    P, μ = positive_projection(pca; f), mean(pca)
    t = P' * (m .- μ)
    t, P, μ
end

## Consistency index

consistency_index(r) = log(sum(>(0), r) / sum(<(0), r))