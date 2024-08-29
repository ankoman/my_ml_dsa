 # export PYTHONPATH="dilithium-py/src:$PYTHONPATH"
from __future__ import annotations
from dilithium_py.ml_dsa import ML_DSA_44
from test_vectors import *
import hashlib

### Testvector generation
print("### Test vector generation ###")
pk, sk = ML_DSA_44.keygen()
# print(pk)
# msg = b"Your message signed by ML_DSA"
# sig = ML_DSA_44.sign(sk, msg)
# assert ML_DSA_44.verify(pk, msg, sig)

# # Verification will fail with the wrong msg or pk
# assert not ML_DSA_44.verify(pk, b"", sig)
# pk_new, sk_new = ML_DSA_44.keygen()
# assert not ML_DSA_44.verify(pk_new, msg, sig)

def i2b(val: int) -> bytearray:
    return val.to_bytes((val.bit_length()+7)//8, 'big')

def b2i(val: bytearray) -> str:
    return int.from_bytes(val, 'big')

def printb(ba: bytearray):
    print(f'{b2i(ba):x}')

def hash_H(din: bytearray, n_bytes: int) -> bytearray:
    s = hashlib.shake_256()
    s.update(din)
    return s.digest(n_bytes)

class polyRing:
    def __init__(self):
        self.q = 8380417
        self.n = 256

    @classmethod
    def rejNTTPoly(cls, rho: bytearray) -> Tq:
        G = hashlib.shake_128()
        G.update(rho)
        G.digest(128)



class my_ml_dsa:
    def __init__(self):
        ### ML_DSA_44
        self.q = 8380417
        self.d = 13  # number of bits dropped from t
        self.tau = 39  # number of ±1 in c
        self.gamma_1 = 131072  # coefficient range of y: 2^17
        self.gamma_2 =  95232  # low order rounding range: (q-1)/88
        self.k = 4  # Dimensions of A = (k, l)
        self.l = 4  # Dimensions of A = (k, l)
        self.eta = 2  # Private key range
        self.omega = 80  # Max number of ones in hint
        self.c_tilde_bytes = 32

    def expandA(self, rho: bytearray):
        for r in range(self.k):
            for s in range(self.l):
                rho_p = rho + s.to_bytes(1, 'little') + r.to_bytes(1, 'little')
                polyRing.rejNTTPoly(rho_p)

        return 0
    
    def _keygen_internal(self, xi: int):
        hash_in = i2b(xi) + self.k.to_bytes(1, 'little') + self.l.to_bytes(1, 'little')
        seeds = hash_H(hash_in, 128)

        rho, rho_p, K = seeds[:32], seeds[32:96], seeds[96:]
        assert b2i(rho) == tv_rho, f'{rho} != {tv_rho:x}'
        assert b2i(rho_p) == tv_rho_p, f'{rho_p} != {tv_rho_p:x}' 
        assert b2i(K) == tv_K, f'{K} != {tv_K:x}' 

        self.expandA(rho)



        return pk, sk

    def keygen(self):
        xi = tv_xi
        pk, sk = self._keygen_internal(xi)


inst = my_ml_dsa()
inst.keygen()