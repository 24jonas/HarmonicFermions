import sqlite3
import numpy as np

def init_db(db_name='recursions_test.db'):
    conn = sqlite3.connect(db_name, timeout=60.0)
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS log_Q_states
                 (n INTEGER, b REAL, m_start INTEGER, L_prev BLOB, L BLOB, log_poch REAL,
                 PRIMARY KEY (n, b))''')
    conn.commit()
    conn.close()

def save_state(n, b, m_start, L_prev, L, log_poch, db_name='recursions_test.db'):
    conn = sqlite3.connect(db_name, timeout=60.0)
    c = conn.cursor()
    c.execute('''INSERT OR REPLACE INTO log_Q_states (n, b, m_start, L_prev, L, log_poch)
                 VALUES (?, ?, ?, ?, ?, ?)''', 
              (n, b, m_start, L_prev.tobytes(), L.tobytes(), log_poch))
    conn.commit()
    conn.close()

def load_state(n, b, db_name='recursions_test.db'):
    conn = sqlite3.connect(db_name, timeout=60.0)
    c = conn.cursor()
    c.execute('''SELECT m_start, L_prev, L, log_poch FROM log_Q_states WHERE b=? AND m_start <= ? ORDER BY m_start DESC LIMIT 1''', (b, n + 2))
    row = c.fetchone()
    conn.close()
    if row is None:
        return None
    m_start = row[0]
    L_prev_old = np.frombuffer(row[1], dtype=np.float64)
    L_old = np.frombuffer(row[2], dtype=np.float64)
    log_poch = row[3]
    
    L_prev = np.zeros(n + 2, dtype=np.float64)
    L = np.zeros(n + 2, dtype=np.float64)
    
    L_prev[:len(L_prev_old)] = L_prev_old
    L[:len(L_old)] = L_old
    
    return m_start, L_prev, L, log_poch

init_db()
# Save a state for n=5
old_n = 5
old_m_start = old_n + 2 # 7
L_prev_old = np.array([1., 2., 3., 4., 5., 6., 7.])
L_old = np.array([7., 6., 5., 4., 3., 2., 1.])
save_state(old_n, 0.5, old_m_start, L_prev_old, L_old, 10.0)

# Load state for n=7
state = load_state(7, 0.5)
if state is not None:
    m_start, L_prev, L, log_poch = state
    print(f"Loaded m_start={m_start}")
    print(f"L_prev shape: {L_prev.shape}, values: {L_prev}")
    print(f"L shape: {L.shape}, values: {L}")

