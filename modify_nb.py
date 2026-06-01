import json

file_path = '/home/jonas/Documents/GitHub/HarmonicFermions/BaxterMeth-opt_freeTaus.ipynb'

with open(file_path, 'r') as f:
    nb = json.load(f)

for cell in nb.get('cells', []):
    if cell.get('cell_type') == 'code':
        source = cell.get('source', [])
        
        new_source = []
        skip = False
        for line in source:
            if 'def plot_fd_vs_tau_dual(' in line:
                line = line.replace('tau_start, tau_end, tau_step', 'tau_list, tau_step')
            
            if 'if tau_end <= tau_start:' in line:
                skip = True
                
            if skip:
                if 'taus_mid = taus_left + 0.5 * tau_step' in line:
                    skip = False
                    new_source.append('    taus_mid = np.array(tau_list)\n')
                continue
                
            if 'tau_start=0.001,' in line:
                line = line.replace('tau_start=0.001', 'tau_list=[1, 0.5, 0.25, 0.125]')
            if 'tau_end=0.25,' in line:
                continue
                
            new_source.append(line)
            
        cell['source'] = new_source

with open(file_path, 'w') as f:
    json.dump(nb, f, indent=1)
    f.write('\n')

print("Notebook modified successfully.")
