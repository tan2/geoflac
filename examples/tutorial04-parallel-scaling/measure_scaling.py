import os
import subprocess
import time
import sys

def main():
    # Overwrite plastic.inp with the 50x150 scaled version
    inp_path = "plastic.inp"
    with open(inp_path, "r") as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        if line.strip().startswith("10,30"):
            new_lines.append("50,150\n")
        elif " 11 " in line:
            new_lines.append(line.replace(" 11 ", " 51 "))
        elif ", 11" in line:
            new_lines.append(line.replace(", 11", ", 51"))
        else:
            new_lines.append(line)

    with open(inp_path, "w") as f:
        f.writelines(new_lines)

    # Thread counts to test
    threads_to_test = [1, 2, 4, 8, 12, 16]
    results = []

    print("==========================================================")
    # Get nproc
    nproc = len(os.sched_getaffinity(0))
    print(f"Starting parallel scaling measurements on {nproc}-core CPU...")
    print("==========================================================")

    # Run for 1 thread first to get baseline
    baseline_time = 0.0

    for num_threads in threads_to_test:
        print(f"Running simulation with {num_threads} threads...", end="", flush=True)
        
        # Clean output files
        for file in os.listdir("."):
            if file.endswith('.0') or file.endswith('.rs') or file.startswith('_contents') or file.startswith('_markers') or file in ['pisos.rs', 'time.rs', 'vbc.s', 'output.asc', 'sys.msg']:
                try:
                    os.remove(file)
                except:
                    pass
                    
        # Set environment and run
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(num_threads)
        
        t0 = time.time()
        res = subprocess.run(["../../src/flac", "plastic.inp"], capture_output=True, env=env)
        t1 = time.time()
        
        elapsed = t1 - t0
        print(f" Done. Time: {elapsed:.2f}s")
        
        if num_threads == 1:
            baseline_time = elapsed
            
        speedup = baseline_time / elapsed
        efficiency = (speedup / num_threads) * 100.0
        
        results.append((num_threads, elapsed, speedup, efficiency))

    # Print results in markdown format
    print("\nScaling Results Table:")
    print("| Threads | Time (s) | Speedup | Efficiency |")
    print("|---------|----------|---------|------------|")
    for r in results:
        print(f"| {r[0]:7d} | {r[1]:8.2f} | {r[2]:7.2f}x | {r[3]:9.1f}% |")

    # Write results to scaling_results.md
    with open("scaling_results.md", "w") as f:
        f.write("# Parallel Scaling Results\n\n")
        f.write(f"Measured on a {nproc}-core CPU using a $50 \\times 150$ grid.\n\n")
        f.write("| Threads | Time (s) | Speedup | Efficiency |\n")
        f.write("|:-------:|:--------:|:-------:|:----------:|\n")
        for r in results:
            f.write(f"| {r[0]} | {r[1]:.2f} | {r[2]:.2f}x | {r[3]:.1f}% |\n")

if __name__ == "__main__":
    main()
