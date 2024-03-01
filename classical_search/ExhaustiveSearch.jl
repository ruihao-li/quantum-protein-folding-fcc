#= 
This program aims to do an exhaustive search of all non-backtracking conformers of a given amino acid sequence of length N and compute the corresponding energy based on a given energy function (Hamiltonian). To implement the non-backtracking constraint, we use the relative turn encoding {0, 1, 2} and then convert it to the absolute turn encoding {0, 1, 2, 3} and then to the full qubit representation {00, 01, 10, 11}. The first two absolute turns are fixed to be [1, 0] and the third turn can only be 1 or 3. The relative turns start from the third turn. So the length of the relative turns is N - 3, where N is the length of the amino acid sequence.
=#
using Base.Threads
using ThreadSafeDicts

# Using all available threads
nthreads()

# Add path
push!(LOAD_PATH, "~/protein-folding-qc/classical_search/")

# println(@__DIR__)

function parse_Hamiltonian(txt_file::String)
    #=
        This function parses the Hamiltonian from a text file and returns the (masked) operators as bitvectors and the corresponding coefficients as floats. The text file should be in the following format:
            + 1.2... * IIIIZI...
            + 0.5... * IIZZI...
            - 0.3... * IZIZI...
            ...
    =#
    coeffs = Vector{Float64}()
    op_arrs = Vector{BitVector}()
    f = open(joinpath(@__DIR__, txt_file), "r+") # open the file
    for line in eachline(f)
        coeff_str, op_str = split(line, " * ")
        if coeff_str[1] == '+'
            coeff = parse(Float64, coeff_str[3:end])
        elseif coeff_str[1] == '-'
            coeff = -parse(Float64, coeff_str[3:end])
        end
        # convert the operator string to a BitArray where 0 = I and 1 = Z
        op_arr = BitArray([c == 'Z' for c in op_str])
        push!(coeffs, coeff)
        push!(op_arrs, op_arr)
    end
    close(f)
    return coeffs, op_arrs
end

@inline _bitwise_AND(a::BitVector, b::BitVector) = a .& b

function compute_energy(coeffs::Vector{Float64}, op_arrs::Vector{BitVector}, config::BitVector)
    #=
    Compute the energy of a given qubit configuration based on the Hamiltonian that consists of only Pauli-Z and identity operators. Bitarrays are used to speed up the computation.
    =#
    energy = 0.0
    @fastmath @simd for i in eachindex(coeffs)
        # bitwise AND between the operator and the configuration (only Z acting on 1 induces a sign change) and get the Hamming weight of sign change
        @inbounds hamming_vec = _bitwise_AND(op_arrs[i], config)
        # determine if counts is odd or even (odd numbers' least significant bit is 1)
        counts = sum(hamming_vec) & 1
        # compute the energy
        energy += counts == 0 ? coeffs[i] : -coeffs[i]
    end
    return energy
end

@inline generate_all_rel_turns(N::Int) = reshape(collect(Iterators.product([0, 2], Iterators.repeated(0:2, N - 4)...)), (2 * 3^(N - 4)))

function rel_turn_to_config_qubits(rel_turns::Tuple{Vararg{Int}})
    #=
    Convert relative turns {0, 1, 2} to absolute turns {0, 1, 2, 3} and then to full qubit representation {00, 01, 10, 11}.
    Note that the first two absolute turns are fixed to be [1, 0] and the third absolute turn can only be 1 or 3. The relative turns start from the third turn. So the length of the relative turns is N - 2, where N is the length of the amino acid sequence. 
    current turn (abs) = (previous turn (abs) + relative turn (rel) + 1) mod 4
    =#
    abs_turns = [1, 0]
    if rel_turns[1] == 1
        error("The first relative turn cannot be 1.")
    end
    for rt in rel_turns
        push!(abs_turns, (abs_turns[end] + rt + 1) % 4)
    end
    # convert to qubit representation as a BitArray (0 -> 00, 1 -> 01, 2 -> 10, 3 -> 11)
    qubit_rep = BitArray(vcat((reverse.(digits.(abs_turns, base=2, pad=2))...)))
    return qubit_rep
end

@inline turn_ind_func_0(i::Int, full_config_seq::BitVector) = ~full_config_seq[2*i-1] & ~full_config_seq[2*i]
@inline turn_ind_func_1(i::Int, full_config_seq::BitVector) = ~full_config_seq[2*i-1] & full_config_seq[2*i]
@inline turn_ind_func_2(i::Int, full_config_seq::BitVector) = full_config_seq[2*i-1] & ~full_config_seq[2*i]
@inline turn_ind_func_3(i::Int, full_config_seq::BitVector) = full_config_seq[2*i-1] & full_config_seq[2*i]

function compute_distance(i::Int, j::Int, full_config_seq::BitVector)
    #=
    This function computes the distance between the i-th and j-th turns of the amino acid sequence.
    =#
    if i > length(full_config_seq) / 2 + 1 || j > length(full_config_seq) / 2 + 1
        error("The turn index is out of range.")
    end
    if i == j
        error("The two turns cannot be the same.")
    end
    if i > j
        i, j = j, i
    end
    dn1 = sum((-1)^k * turn_ind_func_0(k, full_config_seq) for k in i:j-1)
    dn2 = sum((-1)^k * turn_ind_func_1(k, full_config_seq) for k in i:j-1)
    dn3 = sum((-1)^k * turn_ind_func_2(k, full_config_seq) for k in i:j-1)
    dn4 = sum((-1)^k * turn_ind_func_3(k, full_config_seq) for k in i:j-1)
    distance = dn1^2 + dn2^2 + dn3^2 + dn4^2
    return distance
end

function generate_final_qubit_seq(full_config_seq::BitVector)
    #=
    This function generates the final qubit sequence (configuration + interaction) based on the full configuration sequence.
    =#
    full_seq = full_config_seq
    len_peptide = trunc(Int, length(full_config_seq) / 2) + 1
    for i in 1:len_peptide-4
        for j in i+5:len_peptide
            if (j - i) % 2 == 0
                continue
            end
            distance = compute_distance(i, j, full_seq)
            append!(full_seq, distance == 1 ? [true] : [false])
        end
    end
    # reverse the sequence and remove the last four qubits
    full_seq = reverse!(full_seq)[1:end-4]
    # remove the qubit to the left of the last qubit
    deleteat!(full_seq, length(full_seq) - 1)
    return full_seq
end

# Full program
coeffs, op_arrs = parse_Hamiltonian("hamiltonian_7AA_zika.txt");
N = 7; # length of the amino acid sequence
# generate all possible 3-turn encoding {0, 1, 2} of length N - 3
all_rel_turns = generate_all_rel_turns(N)

results = ThreadSafeDict{BitVector,Float64}()
@time @threads for rel_turns in all_rel_turns
    # convert relative turns to absolute turns and then to qubit representation
    full_config_seq = rel_turn_to_config_qubits(rel_turns)
    # generate the final qubit sequence
    final_qubit_seq = generate_final_qubit_seq(full_config_seq)
    # compute the energy
    energy = compute_energy(coeffs, op_arrs, final_qubit_seq)
    # store the result
    results[final_qubit_seq] = energy
end

# write the results to a text file
f = open(joinpath(@__DIR__, "res/", "results_7AA_zika.txt"), "w+")
for (config, energy) in sort(collect(results), by=x -> x[2])[1:50]
    config_str = join(string.(Int.(config)), "")
    write(f, "$config_str $energy\n")
end
close(f)
