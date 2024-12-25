import os
import time
from concurrent.futures import ThreadPoolExecutor

# Node structure
class Node:
    def __init__(self, id, label):
        self.id = id
        self.label = label

# Edge structure
class Edge:
    def __init__(self, from_node, to_node, cigar, variation_type):
        self.from_node = from_node
        self.to_node = to_node
        self.cigar = cigar
        self.variation_type = variation_type

# Global variables for nodes and edges
V = []  # List of nodes
E = []  # List of edges

# Function to read sequences from FASTA file
def read_fasta(filename):
    try:
        with open(filename, "r") as file:
            sequence = ""
            for line in file:
                if line.startswith(">") or not line.strip():
                    continue  # Skip headers and empty lines
                sequence += line.strip()
            return sequence
    except FileNotFoundError:
        print(f"Error: Unable to open FASTA file '{filename}'.")
        return ""

# Function to process a single line of the PAF file
def process_paf_line(line):
    parts = line.strip().split()
    if len(parts) < 8:
        return None, None
    query_start = int(parts[1])
    query_end = int(parts[2])
    variation = parts[7]  # Assuming variation type is in the 8th column
    return [query_start, query_end], variation

# Function to read alignments from PAF file using multithreading
def read_paf(filename):
    breakpoints = []
    variations = []
    try:
        with open(filename, "r") as file:
            lines = file.readlines()
        with ThreadPoolExecutor() as executor:
            results = executor.map(process_paf_line, lines)
            for bp, var in results:
                if bp:
                    breakpoints.extend(bp)
                if var:
                    variations.append(var)
    except FileNotFoundError:
        print(f"Error: Unable to open PAF file '{filename}'.")
    return breakpoints, variations

# Add nodes based on sequence and breakpoints
def process_sequence(sequence, breakpoints):
    global V
    node_id = 1
    prev = 0
    for bp in breakpoints:
        if bp > prev:
            V.append(Node(node_id, sequence[prev:bp]))
            node_id += 1
            prev = bp
    # Capture remaining sequence if any
    if prev < len(sequence):
        V.append(Node(node_id, sequence[prev:]))

# Create edges based on variation information
def create_edges(variations):
    global E
    for i in range(1, len(V)):
        variation_type = variations[i - 1] if i - 1 < len(variations) else "M"
        E.append(Edge(V[i - 1].id, V[i].id, "0M", variation_type))

# Output graph in GFA format
def output_gfa(filename):
    try:
        with open(filename, "w") as file:
            file.write("H\tVN:Z:1.0\n")  # Header
            # Write nodes
            for node in V:
                file.write(f"S\t{node.id}\t{node.label}\n")
            # Write edges
            for edge in E:
                file.write(f"L\t{edge.from_node}\t+\t{edge.to_node}\t+\t{edge.cigar}\t{edge.variation_type}\n")
        print(f"Graph output written to '{filename}'")
    except Exception as e:
        print(f"Error: Unable to write to GFA file '{filename}'. {e}")

# Main function
def main():
    fasta_file = "merged.fasta"
    paf_file = "combined.paf"
    gfa_output = "complex_graph_output.gfa"

    start_time = time.time()

    # Read FASTA file
    sequence = read_fasta(fasta_file)
    if not sequence:
        return

    # Read PAF file
    breakpoints, variations = read_paf(paf_file)

    # Remove duplicate breakpoints and sort
    breakpoints = sorted(set(breakpoints))

    # Process sequence and create nodes
    process_sequence(sequence, breakpoints)

    # Create edges based on variations
    create_edges(variations)

    # Output GFA file
    output_gfa(gfa_output)

    end_time = time.time()
    print(f"Time taken to generate GFA file: {end_time - start_time:.6f} seconds")

if __name__ == "__main__":
    main()
