#!/usr/bin/env python3

import subprocess
import time
import sys
import os
import argparse
import multiprocessing
import re

pattern = re.compile(r'(?m)^(?:\s*[A-Z][A-Z(*)\d\s]*\n)+\s*Cost\s+(\d+)')

def parse_circuit_output(output):
    """Parse the program output to extract the optimal circuit and its cost using regex."""
    # Pattern to match a group of lines starting with capital letters
    # followed by a Cost line

    circuits = []
    costs = []

    # Find all matches in the text
    for match in pattern.finditer(output):
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
        return None, None

    # Find the circuit with lowest cost
    return circuits[0], costs[0]

def compare_circuits(circuit1, circuit2):
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
        # First run to capture output
        result = subprocess.run([program_path] + args, capture_output=True, text=True)

        # Check if program executed successfully
        if result.returncode != 0:
            print(f"Error running {program_path}:")
            print(result.stderr)
            return None, None, None, None

        # Parse the circuit from the output
        circuit, cost = parse_circuit_output(result.stdout)
        if circuit is None:
            print(f"\nDebug - Failed to parse output from {program_path}")
            return None, None, None, None

        # Now run multiple times to get average runtime
        times = []
        for _ in range(num_times):
            start_time = time.perf_counter()
            subprocess.run([program_path] + args, capture_output=True)
            times.append(time.perf_counter() - start_time)

        avg_time = sum(times) / num_times
        return avg_time, circuit, cost, times
    except Exception as e:
        print(f"Error running {program_path}: {str(e)}")
        return None, None, None, None

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
    print(f"Running {args.times} time{'s' if args.times > 1 else ''} for averaging")
    old_time, old_circuit, old_cost, old_times = run_program(old_program, program_args, args.times)

    print(f"\nRunning new version ({new_program})...")
    print(f"Arguments: {' '.join(program_args)}")
    print(f"Running {args.times} time{'s' if args.times > 1 else ''} for averaging")
    new_time, new_circuit, new_cost, new_times = run_program(new_program, program_args, args.times)

    if old_time is not None and new_time is not None:
        print("\nResults:")
        print(f"Old version: {old_time:.3f} seconds (average of {args.times} runs)")
        print(f"Raw times: {', '.join(f'{t:.3f}s' for t in old_times)}")
        print(f"New version: {new_time:.3f} seconds (average of {args.times} runs)")
        print(f"Raw times: {', '.join(f'{t:.3f}s' for t in new_times)}")
        print(f"Speedup: {old_time/new_time:.2f}x")

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