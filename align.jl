#!/usr/local/bin/julia


function importSeqFile(fn::String)
    f = open(fn)
    raw = readstring(f)
    seqdata = Vector{Char}()

    for ch in raw
        push!(seqdata, ch)
    end

    return seqdata
end

function N(mat, row, col)
    if (row == 1 || col == 1)
        println("error: nucleotide coordinates cannot be 1.")
        return
    end

    adjRow = row - 1

    return mat[adjRow, col] - 2
end

function W(mat, row, col)
    if (row == 1 || col == 1)
        println("error: nucleotide coordinates cannot be 1.")
        return
    end

    adjCol = col - 1

    return mat[row, adjCol] - 2
end

function matchnuc(s, t, row, col)
    if (row == 1 || col == 1)
        println("error: nucleotide coordinates cannot be 1.")
        return
    end

    if (s[row-1] == t[col-1])
        return 1
    else
        return -1;
    end
end

function NW(mat, s, t, row, col)
    if (row == 1 || col == 1)
        println("error: nucleotide coordinates cannot be 1.")
        return
    end

    adjRow = row - 1
    adjCol = col - 1

    return mat[adjRow, adjCol] + matchnuc(s, t, row, col)
end

function source(idx, row, col)
    if (row == 1 || col == 1)
        println("error: nucleotide coordinates cannot be 1.")
        return
    end

    if idx == 1
        row = row - 1
    elseif idx == 2
        row = row - 1
        col = col - 1
    elseif idx == 3
        col = col - 1;
    else
        println("error: invalid input.")
    end

    return (row, col)
end

function SmithWaterman(sim_mat, tup_mat, s, t, row, col)
    if (row == 1 || col == 1)
        println("error: nucleotide coordinates cannot be 1.")
        return
    end

    scores = Vector{Int}()
    push!(scores, N(sim_mat, row, col))
    push!(scores, NW(sim_mat, s, t, row, col))
    push!(scores, W(sim_mat, row, col))
    # @show scores

    top_score = maximum(scores)
    top_index = indmax(scores)

    if top_score < 0
        top_score = 0
    end

    sim_mat[row, col] = top_score

    # tup_mat[row, col] = source(top_index, row, col)
end

function backtrack(mat)
    x_sz = size(mat, 1)
    y_sz = size(mat, 2)

    cur_max = mat[x_sz, y_sz]
    row = 1
    col = 1

    for j in y_sz:-1:2
        for i in x_sz:-1:2
            if mat[i, j] > cur_max
                cur_max = mat[i, j]
                row = i
                col = j
            end
        end
    end

    return (cur_max, [row-1, col-1])
end



# seqFilNam = "ex1_seq.txt"
# unkFilNam = "ex1_unk.txt"

seqFilNam = "ex2_seq.txt"
unkFilNam = "ex2_unk.txt"

# seqFilNam = "HIV-1_db.fasta"
# unkFilNam = "HIV-1_Polymerase.txt"

s = importSeqFile(seqFilNam)
pop!(s)

t = importSeqFile(unkFilNam)
pop!(t)

sim_mat = zeros(Matrix{Int}(length(s)+1,
                            length(t)+1))

for j in 2:size(sim_mat, 2)
    for i in 2:size(sim_mat, 1)
        sim_mat[i, j] = -999
    end
end

#display(sim_mat)

tup_mat = Matrix{Tuple{Int, Int}}(length(s)+1,
                                  length(t)+1)

for j in 1:size(tup_mat, 2)
    for i in 1:size(tup_mat, 1)
        tup_mat[i, j] = (0, 0)
    end
end

# for j in 2:size(sim_mat, 2)
#     for i in 2:size(sim_mat, 1)
#         SmithWaterman(sim_mat, tup_mat, s, t, i, j)
#     end
# end

# display(sim_mat)
println(backtrack(sim_mat))
