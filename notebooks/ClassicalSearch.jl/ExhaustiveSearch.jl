#= 
This program aims to do an exhaustive search of all non-backtracking conformers of a given amino acid sequence of length N and compute the corresponding energy based on a given energy function (Hamiltonian). To implement the non-backtracking constraint, we use the relative turn encoding {0, 1, 2} and then convert it to the absolute turn encoding {0, 1, 2, 3} and then to the full qubit representation {00, 01, 10, 11}. The first two absolute turns are fixed to be [1, 0] and the third turn can only be 1 or 3. The relative turns start from the third turn. So the length of the relative turns is N - 3, where N is the length of the amino acid sequence.
=#
using Base.Threads
using ThreadSafeDicts

# Using all available threads
nthreads()

# Add path
push!(LOAD_PATH, "/Users/ruihaoli/Library/CloudStorage/OneDrive-TheUniversityofSydney (Students)/CCF Work/protein_folding/temp_codes")

println(@__DIR__)

function parse_Hamiltonian(txt_file::String)
    #=
    This function parses the Hamiltonian from a text file and returns a dictionary of the Hamiltonian terms as {operator: coeff}. The text file should be in the following format:
        + 1.2... * IIIIZI...
        + 0.5... * IIZZI...
        - 0.3... * IZIZI...
        ...
    =#
    hamiltonian = Dict()
    f = open(joinpath(@__DIR__, txt_file), "r+") # open the file
    for line in eachline(f)
        coeff_str, op_str = split(line, " * ")
        if coeff_str[1] == '+'
            coeff = parse(Float64, coeff_str[3:end])
        elseif coeff_str[1] == '-'
            coeff = -parse(Float64, coeff_str[3:end])
        end
        hamiltonian[op_str] = coeff
    end
    close(f)
    return hamiltonian
end

# parse_Hamiltonian("hamiltonian_zika.txt")

function compute_energy(hamiltonian::Dict, config::String)
    #=
    Compute the energy of a given qubit configuration based on the Hamiltonian that consists of only Pauli-Z and identity operators.
    =#
    energy = 0.0
    for (op, coeff) in hamiltonian
        # compute the energy of each term
        z_pos = findall(x -> x == 'Z', op)
        # count how many 1s in the corresponding qubit positions
        term_energy = length(z_pos) == 0 ? 1.0 : (-1.0)^sum(parse(Int, config[i]) for i in z_pos)
        energy += coeff * term_energy
    end
    return energy
end

# compute_energy(parse_Hamiltonian("hamiltonian_angiotensin.txt"), "1000100001110001110011")
# compute_energy(parse_Hamiltonian("hamiltonian_zika.txt"), "011000101")

function generate_all_rel_turns(N::Int)
    #=
    This function generates all possible 3-turn (relative turns) encoding {0, 1, 2} of length N - 3, where N is the number of amino acids. Note that the first turn can only be 0 or 2.
    =#
    all_rel_turns = reshape(collect(Iterators.product([0, 2], Iterators.repeated(0:2, N - 4)...)), (2 * 3^(N - 4)))
    # convert to strings ['000...', '002...', ...]
    all_rel_turns = [string(turn...) for turn in all_rel_turns]
    return all_rel_turns
end

# generate_all_rel_turns(7)

function rel_turn_to_config_qubits(rel_turns::String)
    #=
    Convert relative turns {0, 1, 2} to absolute turns {0, 1, 2, 3} and then to full qubit representation {00, 01, 10, 11}.
    Note that the first two absolute turns are fixed to be [1, 0] and the third absolute turn can only be 1 or 3. The relative turns start from the third turn. So the length of the relative turns is N - 2, where N is the length of the amino acid sequence. 
    current turn (abs) = (previous turn (abs) + relative turn (rel) + 1) mod 4
    =#
    # abs_turns = "10"
    qubit_rep = "0100"
    if rel_turns[1] == '1'
        error("The first relative turn cannot be 1.")
    end
    for rt in rel_turns
        rt = parse(Int, rt)
        this_turn = (parse(Int, qubit_rep[end-1:end], base=2) + rt + 1) % 4
        # abs_turns *= string(this_turn)
        qubit_rep *= join(reverse(digits(this_turn, base=2, pad=2)))
    end
    return qubit_rep
end

# rel_turn_to_config_qubits("2022")

function turn_indicator_funcs(i::Int, full_config_seq::String)
    #=
    This function returns the indicator functions for the i-th turn of the full configuration sequence. 
    =#
    if i > length(full_config_seq) / 2
        error("The turn index cannot be larger than the half of the length of the full configuration sequence.")
    end
    f0 = (1 - parse(Int, full_config_seq[2*i-1])) * (1 - parse(Int, full_config_seq[2*i]))
    f1 = parse(Int, full_config_seq[2*i]) * (parse(Int, full_config_seq[2*i]) - parse(Int, full_config_seq[2*i-1]))
    f2 = parse(Int, full_config_seq[2*i-1]) * (parse(Int, full_config_seq[2*i-1]) - parse(Int, full_config_seq[2*i]))
    f3 = parse(Int, full_config_seq[2*i-1]) * parse(Int, full_config_seq[2*i])
    return f0, f1, f2, f3
end

# turn_indicator_funcs(2, "010011010001")

function compute_distance(i::Int, j::Int, full_config_seq::String)
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
    dn1 = sum((-1)^k * turn_indicator_funcs(k, full_config_seq)[1] for k in i:j-1)
    dn2 = sum((-1)^k * turn_indicator_funcs(k, full_config_seq)[2] for k in i:j-1)
    dn3 = sum((-1)^k * turn_indicator_funcs(k, full_config_seq)[3] for k in i:j-1)
    dn4 = sum((-1)^k * turn_indicator_funcs(k, full_config_seq)[4] for k in i:j-1)
    distance = dn1^2 + dn2^2 + dn3^2 + dn4^2
    return distance
end

# compute_distance(2, 5, "010011010001")

function generate_final_qubit_seq(full_config_seq::String)
    #=
    This function generates the final qubit sequence (configuration + interaction) based on the full configuration sequence.
    =#
    full_seq = full_config_seq
    len_peptide = trunc(Int, length(full_seq) / 2) + 1
    for i in 1:len_peptide-4
        for j in i+5:len_peptide
            if (j - i) % 2 == 0
                continue
            end
            distance = compute_distance(i, j, full_seq)
            full_config_seq *= distance == 1 ? "1" : "0"
        end
    end
    # reverse the sequence and remove the last four qubits
    final_seq = reverse(full_config_seq)[1:end-4]
    # remove the qubit to the left of the last qubit
    final_seq = final_seq[1:end-2] * final_seq[end]
    return final_seq
end

# generate_final_qubit_seq(rel_turn_to_config_qubits("2022")[2])
# compute_energy(hamiltonian, generate_final_qubit_seq("010001111001001001"))

# Full program
hamiltonian = parse_Hamiltonian("hamiltonian_13AA_fake.txt");
N = 13; # length of the amino acid sequence
# generate all possible 3-turn encoding {0, 1, 2} of length N - 3
all_rel_turns = generate_all_rel_turns(N)

results = ThreadSafeDict{String,Float64}()
@time @threads for rel_turns in all_rel_turns
    # convert relative turns to absolute turns and then to qubit representation
    full_config_seq = rel_turn_to_config_qubits(rel_turns)
    # generate the final qubit sequence
    final_qubit_seq = generate_final_qubit_seq(full_config_seq)
    # compute the energy
    energy = compute_energy(hamiltonian, final_qubit_seq)
    # write {config: energy} to a dictionary in a threadsafe way
    results[final_qubit_seq] = energy
end

# sort the results by energy and list the 10 lowest energy conformers
sort(collect(results), by=x -> x[2])[1:10]

