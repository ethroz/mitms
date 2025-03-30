#!/usr/bin/env python3

import subprocess
import time
import sys
import os
import argparse
import multiprocessing
import re

circuit_cost_pattern = re.compile(r'(?m)^(?:\s*[A-Z][A-Z(*)\d\s]*\n)+\s*Cost\s+(\d+)')
time_pattern = re.compile(r'Time:\s*(\d+\.\d+)\s*s')
search_pattern = re.compile(r'# searches so far:\s*(\d+)')

def parse_timing_and_searches(output):
    """Parse timing and search count information from the output."""    
    generation_times = []
    searching_times = []
    searches = []
    
    # Find all time entries
    for i, match in enumerate(time_pattern.finditer(output)):
        match_time = float(match.group(1))
        if i == 0 or i % 3 == 2:
            generation_times.append(match_time)
        else:
            searching_times.append(match_time)
    
    # Find all search count entries
    for match in search_pattern.finditer(output):
        searches.append(int(match.group(1)))
    
    return generation_times, searching_times, searches

def parse_circuit_output(output):
    """Parse the program output to extract the optimal circuit and its cost using regex."""
    circuits = []
    costs = []
    generation_times, searching_times, searches = parse_timing_and_searches(output)

    # Find all matches in the text
    for match in circuit_cost_pattern.finditer(output):
        # Get the full matched text
        circuit_text = match.group(0)

        # Split into lines and filter out empty lines
        lines = [line.strip() for line in circuit_text.split('\n') if line.strip()]

        # Last line contains the cost, all others are circuit lines
        circuit_lines = lines[:-1]  # All lines except the last one
        cost = int(lines[-1].replace('Cost', '').strip())

        # Join the circuit lines with newlines
        circuit = '\n'.join(circuit_lines)
        circuits.append(circuit)
        costs.append(cost)

    if not circuits:
        print("\nDebug - No circuits found in output. Raw output:")
        print(output)
        return None, None, None, None

    # Find the circuit with lowest cost
    return circuits[0], costs[0], generation_times, searching_times, searches

def compare_circuits(circuit1: str, circuit2: str):
    """Compare two circuits and return detailed differences."""
    if circuit1 == circuit2:
        return None

    # Split circuits into lines for comparison
    lines1 = circuit1.split('\n')
    lines2 = circuit2.split('\n')

    # Find the first line where they differ
    for i, (line1, line2) in enumerate(zip(lines1, lines2)):
        if line1 != line2:
            return f"First difference at line {i+1}:\nOld: {line1}\nNew: {line2}"

    # If one circuit is longer than the other
    if len(lines1) != len(lines2):
        return f"Different number of lines:\nOld: {len(lines1)} lines\nNew: {len(lines2)} lines"

    return "Unknown difference"

def run_program(program_path, args, num_times=1):
    """Run a program with given arguments and return its execution time and circuit output."""
    try:
        total_times = []
        result = None
        all_generation_times = []
        all_searching_times = []
        all_searches = []

        # Run the program num_times
        for i in range(num_times):
            start_time = time.perf_counter()
            current_result = subprocess.run([program_path] + args, capture_output=True, text=True)
            end_time = time.perf_counter()
            total_times.append(end_time - start_time)
            
            # Check if this run failed
            if current_result.returncode != 0:
                print(f"Error running {program_path} (run {i+1}/{num_times}):")
                print(current_result.stderr)
                return None, None, None, None, None, None, None
            
            # Store the first run's result for circuit parsing
            if i == 0:
                result = current_result
            
            # Parse timing and search information from each run
            gen_times, search_times, searches = parse_timing_and_searches(current_result.stdout)
            all_generation_times.append(gen_times)
            all_searching_times.append(search_times)
            all_searches.append(searches)

        # Parse the circuit from the first run's output
        circuit, cost, _, _, _ = parse_circuit_output(result.stdout)
        if circuit is None:
            print(f"\nDebug - Failed to parse output from {program_path}")
            return None, None, None, None, None, None, None

        # Average the timing and search information across all runs
        avg_generation_times = []
        avg_searching_times = []
        avg_searches = []
        
        # Average generation times
        for i in range(len(all_generation_times[0])):
            avg_generation_times.append(sum(run[i] for run in all_generation_times) / num_times)
            
        # Average searching times
        for i in range(len(all_searching_times[0])):
            avg_searching_times.append(sum(run[i] for run in all_searching_times) / num_times)
            
        # Average search counts
        for i in range(len(all_searches[0])):
            avg_searches.append(round(sum(run[i] for run in all_searches) / num_times))

        avg_time = sum(total_times) / num_times
        return avg_time, circuit, cost, total_times, avg_generation_times, avg_searching_times, avg_searches
    except Exception as e:
        print(f"Error running {program_path}: {str(e)}")
        return None, None, None, None, None, None, None

def main():
    parser = argparse.ArgumentParser(description='Compare runtime of two MITMS program versions')
    parser.add_argument('circuit', help='The quantum circuit to synthesize')
    parser.add_argument('--threads', type=int,
                       default=multiprocessing.cpu_count(),
                       help='Number of threads to use (default: number of CPU cores)')
    parser.add_argument('--times', type=int,
                       default=1,
                       help='Number of times to run each program for averaging (default: 1)')
    args = parser.parse_args()

    # Get the platform-specific paths
    old_program = "./build/mitms_original"
    new_program = "./build/mitms_parallel"

    # Check if both programs exist
    if not os.path.exists(old_program):
        print(f"Error: Old program not found at {old_program}")
        sys.exit(1)
    if not os.path.exists(new_program):
        print(f"Error: New program not found at {new_program}")
        sys.exit(1)
    
    # Check the arguments
    if args.times < 1:
        print("Error: Number of times to run each program must be at least 1")
        sys.exit(1)
    if args.threads < 1:
        print("Error: Number of threads must be at least 1")
        sys.exit(1)

    # Prepare arguments
    program_args = [
        "-no-serialize",
        "-no-equiv-checks",
        "-threads", str(args.threads),
        "-early-stop",
        args.circuit
    ]

    print(f"Running old version ({old_program})...")
    print(f"Arguments: {' '.join(program_args)}")
    print(f"Running {args.times} time{'s for averaging' if args.times > 1 else ''}")
    old_time, old_circuit, old_cost, old_times, old_generation_times, \
        old_searching_times, old_searches = run_program(old_program, program_args, args.times)

    print(f"\nRunning new version ({new_program})...")
    print(f"Arguments: {' '.join(program_args)}")
    print(f"Running {args.times} time{'s for averaging' if args.times > 1 else ''}")
    new_time, new_circuit, new_cost, new_times, new_generation_times, \
        new_searching_times, new_searches = run_program(new_program, program_args, args.times)

    if old_time is not None and new_time is not None:
        print("\nResults:")
        print(f"Old version: {old_time:.3f}s{(f' (average of {args.times} runs)' if args.times > 1 else '')}")
        print(f"Raw times: {', '.join(f'{t:.3f}s' for t in old_times)}")
        print(f"New version: {new_time:.3f}s{(f' (average of {args.times} runs)' if args.times > 1 else '')}")
        print(f"Raw times: {', '.join(f'{t:.3f}s' for t in new_times)}")
        print(f"Speedup: {old_time/new_time:.2f}x")

        print("\nStep-by-step timing comparison:")
        print(f"Old version generation times: {', '.join(f'{t:.3f}s' for t in old_generation_times)}")
        print(f"New version generation times: {', '.join(f'{t:.3f}s' for t in new_generation_times)}")
        print(f"Old version searching times: {', '.join(f'{t:.3f}s' for t in old_searching_times)}")
        print(f"New version searching times: {', '.join(f'{t:.3f}s' for t in new_searching_times)}")
        
        print("\nSearch count comparison:")
        print(f"Old version searches: {', '.join(f'{s}' for s in old_searches)}")
        print(f"New version searches: {', '.join(f'{s}' for s in new_searches)}")

        print("\nCircuit Comparison:")
        if old_circuit and new_circuit:
            print(f"Old version optimal circuit (cost: {old_cost}):")
            print(old_circuit)
            print(f"\nNew version optimal circuit (cost: {new_cost}):")
            print(new_circuit)

            if old_circuit == new_circuit:
                print("\n✓ Circuits are identical")
            else:
                print("\n✗ Circuits are different!")
                diff = compare_circuits(old_circuit, new_circuit)
                if diff:
                    print("\nDetailed difference:")
                    print(diff)
        else:
            print("Could not parse circuits from output")
    else:
        print("\nComparison failed due to errors. See above for details.")

if __name__ == "__main__":
    main() 