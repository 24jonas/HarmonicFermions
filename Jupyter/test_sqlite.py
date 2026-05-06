import sqlite3
import numpy as np

def init_db(db_name="recursions.db"):
    conn = sqlite3.connect(db_name, timeout=60.0)
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS log_Q_states
                 (n INTEGER, b REAL, m_start INTEGER, L_prev BLOB, L BLOB, log_poch REAL,
                 PRIMARY KEY (n, b))''')
    conn.commit()
    conn.close()

def save_log_Q_state(n, b, m_start, L_prev, L, log_poch, db_name="recursions.db"):
    conn = sqlite3.connect(db_name, timeout=60.0)
    c = conn.cursor()
    c.execute('''INSERT OR REPLACE INTO log_Q_states (n, b, m_start, L_prev, L, log_poch)
                 VALUES (?, ?, ?, ?, ?, ?)''', 
              (n, b, m_start, L_prev.tobytes(), L.tobytes(), log_poch))
    conn.commit()
    conn.close()

def load_log_Q_state(n, b, db_name="recursions.db"):
    conn = sqlite3.connect(db_name, timeout=60.0)
    c = conn.cursor()
    c.execute('''SELECT m_start, L_prev, L, log_poch FROM log_Q_states WHERE n=? AND b=?''', (n, b))
    row = c.fetchone()
    conn.close()
    if row is None:
        return None
    m_start = row[0]
    L_prev = np.frombuffer(row[1], dtype=np.float64).copy()
    L = np.frombuffer(row[2], dtype=np.float64).copy()
    log_poch = row[3]
    return m_start, L_prev, L, log_poch

init_db()
L = np.zeros(10)
save_log_Q_state(10, 0.5, 2, L, L, 0.0)
res = load_log_Q_state(10, 0.5)
print("Loaded:", res[0], type(res[1]), res[1].flags.writeable)
