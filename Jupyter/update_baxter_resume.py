import json

notebook_path = "/home/jonas/Documents/GitHub/HarmonicFermions/Jupyter/BaxterMeth-opt.ipynb"

with open(notebook_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

for i, cell in enumerate(nb['cells']):
    src = "".join(cell['source'])
    if "def load_log_Q_state(" in src:
        # We replace the two load functions
        new_src = src.replace(
            "c.execute('''SELECT m_start, L_prev, L, log_poch FROM log_Q_states WHERE n=? AND b=?''', (n, b))",
            "c.execute('''SELECT m_start, L_prev, L, log_poch FROM log_Q_states WHERE b=? AND m_start <= ? ORDER BY m_start DESC LIMIT 1''', (b, n + 2))"
        )
        
        # Then replace the assignment and padding
        old_assign = """            m_start = row[0]
            L_prev = np.frombuffer(row[1], dtype=np.float64).copy()
            L = np.frombuffer(row[2], dtype=np.float64).copy()
            log_poch = row[3]
            return m_start, L_prev, L, log_poch"""
            
        new_assign = """            m_start = row[0]
            L_prev_old = np.frombuffer(row[1], dtype=np.float64)
            L_old = np.frombuffer(row[2], dtype=np.float64)
            log_poch = row[3]
            L_prev = np.zeros(n + 2, dtype=np.float64)
            L = np.zeros(n + 2, dtype=np.float64)
            L_prev[:len(L_prev_old)] = L_prev_old
            L[:len(L_old)] = L_old
            return m_start, L_prev, L, log_poch"""
        new_src = new_src.replace(old_assign, new_assign)
        
        # Now for load_log_Q_all_state
        new_src = new_src.replace(
            "c.execute('''SELECT m_start, L_prev, L, res, log_poch FROM log_Q_all_states WHERE n=? AND b=?''', (n, b))",
            "c.execute('''SELECT m_start, L_prev, L, res, log_poch FROM log_Q_all_states WHERE b=? AND m_start <= ? ORDER BY m_start DESC LIMIT 1''', (b, n + 2))"
        )
        
        old_assign_all = """            m_start = row[0]
            L_prev = np.frombuffer(row[1], dtype=np.float64).copy()
            L = np.frombuffer(row[2], dtype=np.float64).copy()
            res = np.frombuffer(row[3], dtype=np.float64).copy()
            log_poch = row[4]
            return m_start, L_prev, L, res, log_poch"""
            
        new_assign_all = """            m_start = row[0]
            L_prev_old = np.frombuffer(row[1], dtype=np.float64)
            L_old = np.frombuffer(row[2], dtype=np.float64)
            res_old = np.frombuffer(row[3], dtype=np.float64)
            log_poch = row[4]
            L_prev = np.zeros(n + 2, dtype=np.float64)
            L = np.zeros(n + 2, dtype=np.float64)
            res = np.zeros(n + 1, dtype=np.float64)
            L_prev[:len(L_prev_old)] = L_prev_old
            L[:len(L_old)] = L_old
            res[:len(res_old)] = res_old
            return m_start, L_prev, L, res, log_poch"""
            
        new_src = new_src.replace(old_assign_all, new_assign_all)
        
        # Replace the cell source
        # Note we have to split it back into lines for the JSON
        lines = [line + '\n' for line in new_src.split('\n')]
        # The last split item will be empty string and have a newline appended, we can fix it if needed but Jupyter doesn't care
        lines[-1] = lines[-1].replace('\n', '')
        nb['cells'][i]['source'] = lines

with open(notebook_path, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)

print("BaxterMeth-opt updated successfully to resume from ANY smaller n.")
