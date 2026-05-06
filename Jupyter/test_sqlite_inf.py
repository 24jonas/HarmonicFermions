import sqlite3
import numpy as np

conn = sqlite3.connect(':memory:')
c = conn.cursor()
c.execute('CREATE TABLE test (val REAL)')
c.execute('INSERT INTO test VALUES (?)', (-np.inf,))
c.execute('SELECT * FROM test')
print("Returned:", c.fetchone()[0])
