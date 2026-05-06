import json

notebook_path = "/home/jonas/Documents/GitHub/HarmonicFermions/Jupyter/JonasPolyConv-opt.ipynb"

with open(notebook_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

for i, cell in enumerate(nb['cells']):
    src = "".join(cell['source'])
    if "def load_polyconv_state(" in src:
        # We replace the query
        new_src = src.replace(
            "c.execute('''SELECT k_start, logZ, prev_logZn, small_count, cumulative_capacity FROM polyconv_states WHERE n=? AND b=? AND d=?''', (n, b, d))",
            "c.execute('''SELECT n, k_start, logZ, prev_logZn, small_count, cumulative_capacity FROM polyconv_states WHERE d=? AND b=? AND n >= ? ORDER BY k_start DESC LIMIT 1''', (d, b, n))"
        )
        
        # Then replace the assignment and slicing
        old_assign = """            k_start = row[0]
            logZ = np.frombuffer(row[1], dtype=np.float64).copy()
            prev_logZn = row[2]
            small_count = row[3]
            cumulative_capacity = row[4]
            return k_start, logZ, prev_logZn, small_count, cumulative_capacity"""
            
        new_assign = """            saved_n = row[0]
            k_start = row[1]
            logZ_saved = np.frombuffer(row[2], dtype=np.float64)
            prev_logZn = row[3]
            small_count = row[4]
            cumulative_capacity = row[5]
            logZ = logZ_saved[:n + 1].copy()
            return k_start, logZ, prev_logZn, small_count, cumulative_capacity"""
            
        new_src = new_src.replace(old_assign, new_assign)
        
        lines = [line + '\n' for line in new_src.split('\n')]
        lines[-1] = lines[-1].replace('\n', '')
        nb['cells'][i]['source'] = lines

with open(notebook_path, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)

print("JonasPolyConv-opt updated successfully to resume from ANY larger n.")
