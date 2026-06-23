"""
Baxter Method — Partition function computation via q-series recursion.

Computes log Q_n(b) for n fermions in a 2D harmonic trap using Baxter's
combinatorial identity. Supports single-n (log_Q) and all-n (log_Q_all)
evaluation, with SQLite-based checkpointing for long computations.

References:
    Baxter, R. J. (1982). Exactly Solved Models in Statistical Mechanics.
"""

import math
import numpy as np
import numba
import sqlite3
import time


# ============================================================
# Database persistence (checkpoint / resume long computations)
# ============================================================

def init_db(db_name='recursions.db'):
    """Create the SQLite tables for checkpointing if they don't exist."""
    try:
        conn = sqlite3.connect(db_name, timeout=60.0)
        c = conn.cursor()
        c.execute('''CREATE TABLE IF NOT EXISTS log_Q_states
                     (n INTEGER, b REAL, m_start INTEGER, L_prev BLOB, L BLOB, log_poch REAL,
                     PRIMARY KEY (n, b))''')
        c.execute('''CREATE TABLE IF NOT EXISTS log_Q_all_states
                     (n INTEGER, b REAL, m_start INTEGER, L_prev BLOB, L BLOB, res BLOB, log_poch REAL,
                     PRIMARY KEY (n, b))''')
        conn.commit()
        conn.close()
    except Exception:
        pass


def save_log_Q_state(n, b, m_start, L_prev, L, log_poch, db_name='recursions.db'):
    """Save a single-n computation checkpoint."""
    for _ in range(5):
        try:
            conn = sqlite3.connect(db_name, timeout=60.0)
            c = conn.cursor()
            c.execute('''INSERT OR REPLACE INTO log_Q_states (n, b, m_start, L_prev, L, log_poch)
                         VALUES (?, ?, ?, ?, ?, ?)''',
                      (n, b, m_start, L_prev.tobytes(), L.tobytes(), log_poch))
            conn.commit()
            conn.close()
            break
        except sqlite3.OperationalError:
            time.sleep(2)


def load_log_Q_state(n, b, db_name='recursions.db'):
    """Load a single-n computation checkpoint, or return None."""
    for _ in range(5):
        try:
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
        except sqlite3.OperationalError:
            time.sleep(2)
    return None


def save_log_Q_all_state(n, b, m_start, L_prev, L, res, log_poch, db_name='recursions.db'):
    """Save an all-n computation checkpoint."""
    for _ in range(5):
        try:
            conn = sqlite3.connect(db_name, timeout=60.0)
            c = conn.cursor()
            c.execute('''INSERT OR REPLACE INTO log_Q_all_states (n, b, m_start, L_prev, L, res, log_poch)
                         VALUES (?, ?, ?, ?, ?, ?, ?)''',
                      (n, b, m_start, L_prev.tobytes(), L.tobytes(), res.tobytes(), log_poch))
            conn.commit()
            conn.close()
            break
        except sqlite3.OperationalError:
            time.sleep(2)


def load_log_Q_all_state(n, b, db_name='recursions.db'):
    """Load an all-n computation checkpoint, or return None."""
    for _ in range(5):
        try:
            conn = sqlite3.connect(db_name, timeout=60.0)
            c = conn.cursor()
            c.execute('''SELECT m_start, L_prev, L, res, log_poch FROM log_Q_all_states WHERE b=? AND m_start <= ? ORDER BY m_start DESC LIMIT 1''', (b, n + 2))
            row = c.fetchone()
            conn.close()
            if row is None:
                return None
            m_start = row[0]
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
            return m_start, L_prev, L, res, log_poch
        except sqlite3.OperationalError:
            time.sleep(2)
    return None


# ============================================================
# Numba-JIT computational kernels
# ============================================================

@numba.njit
def logsumexp2(x, y):
    """Log-sum-exp of two values."""
    m = x if x > y else y
    return m + math.log(math.exp(x - m) + math.exp(y - m))


@numba.njit
def logsumexp_array(vals):
    """Log-sum-exp over an array."""
    m = np.max(vals)
    s = 0.0
    for i in range(len(vals)):
        s += math.exp(vals[i] - m)
    return m + math.log(s)


@numba.njit
def log1mexp_neg(x):
    """Compute log(1 - exp(-|x|)) = log(-expm1(x)) for x < 0."""
    return math.log(-math.expm1(x))


@numba.njit
def _log_Q_chunk(m_start, m_end, logb, L_prev, L):
    """Process a chunk of the recursion for single-n log_Q."""
    for m in range(m_start, m_end):
        L[m - 1] = logsumexp_array(L_prev[:m - 1])
        log1mbpow = log1mexp_neg((m - 1) * logb)
        for i in range(m - 1, 0, -1):
            t1 = logb + L[i]
            t2 = (1 - i) * logb + log1mbpow + L_prev[i - 1]
            L[i - 1] = logsumexp2(t1, t2)
        temp = L_prev
        L_prev = L
        L = temp
    return L_prev, L


@numba.njit
def _compute_log_poch(n, logb):
    """Compute log of the q-Pochhammer symbol."""
    log_poch = 0.0
    for k in range(1, n + 1):
        log_poch += log1mexp_neg(k * logb)
    return log_poch


@numba.njit
def _log_Q_all_chunk(m_start, m_end, logb, L_prev, L, res, log_poch):
    """Process a chunk of the recursion, storing results for all n."""
    for m in range(m_start, m_end):
        L[m - 1] = logsumexp_array(L_prev[:m - 1])
        log1mbpow = log1mexp_neg((m - 1) * logb)
        for i in range(m - 1, 0, -1):
            t1 = logb + L[i]
            t2 = (1 - i) * logb + log1mbpow + L_prev[i - 1]
            L[i - 1] = logsumexp2(t1, t2)
        temp = L_prev
        L_prev = L
        L = temp
        num = m - 1
        logH_num = L_prev[num]
        log_poch += log1mexp_neg(num * logb)
        res[num] = (num * (num + 1) / 2.0) * logb - 2.0 * log_poch + logH_num
    return L_prev, L, res, log_poch


# ============================================================
# Main API
# ============================================================

def log_Q(n, b, db_name='recursions.db'):
    """
    Compute log Q_n(b) — the log partition function for n fermions.

    Parameters
    ----------
    n : int
        Number of fermions.
    b : float
        Boltzmann factor, must be in (0, 1).
    db_name : str
        Path to the SQLite database for checkpointing.

    Returns
    -------
    float
        The value log Q_n(b).
    """
    if not (0.0 < b < 1.0):
        raise ValueError("b must lie in (0,1)")
    logb = math.log(b)
    init_db(db_name)
    state = load_log_Q_state(n, b, db_name)
    if state is not None:
        m_start, L_prev, L, log_poch = state
    else:
        m_start = 2
        L_prev = np.zeros(n + 2, dtype=np.float64)
        L = np.zeros(n + 2, dtype=np.float64)
        log_poch = 0.0

    chunk_size = 50000
    while m_start < n + 2:
        m_end = min(m_start + chunk_size, n + 2)
        L_prev, L = _log_Q_chunk(m_start, m_end, logb, L_prev, L)
        m_start = m_end
        save_log_Q_state(n, b, m_start, L_prev, L, log_poch, db_name)

    logH = L_prev[n]
    log_poch = _compute_log_poch(n, logb)
    return (n * (n + 1) // 2) * logb - 2.0 * log_poch + logH


def log_Q_all(n, b, db_name='recursions.db'):
    """
    Compute log Q_k(b) for all k = 0, 1, ..., n.

    Parameters
    ----------
    n : int
        Maximum number of fermions.
    b : float
        Boltzmann factor, must be in (0, 1).
    db_name : str
        Path to the SQLite database for checkpointing.

    Returns
    -------
    np.ndarray
        Array of shape (n+1,) with res[k] = log Q_k(b).
    """
    if not (0.0 < b < 1.0):
        raise ValueError("b must lie in (0,1)")
    logb = math.log(b)
    init_db(db_name)
    state = load_log_Q_all_state(n, b, db_name)
    if state is not None:
        m_start, L_prev, L, res, log_poch = state
    else:
        m_start = 2
        L_prev = np.zeros(n + 2, dtype=np.float64)
        L = np.zeros(n + 2, dtype=np.float64)
        res = np.zeros(n + 1, dtype=np.float64)
        log_poch = 0.0

    chunk_size = 50000
    while m_start < n + 2:
        m_end = min(m_start + chunk_size, n + 2)
        L_prev, L, res, log_poch = _log_Q_all_chunk(m_start, m_end, logb, L_prev, L, res, log_poch)
        m_start = m_end
        save_log_Q_all_state(n, b, m_start, L_prev, L, res, log_poch, db_name)

    return res
