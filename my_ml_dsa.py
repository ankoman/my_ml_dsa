 # export PYTHONPATH="dilithium-py/src:$PYTHONPATH"
from __future__ import annotations
from typing import List
from dilithium_py.ml_dsa import ML_DSA_44
from test_vectors import *
import hashlib, copy

zetas = [0, 4808194, 3765607, 3761513, 5178923, 5496691, 5234739, 5178987, 7778734, 3542485, 2682288, 2129892, 3764867, 7375178, 557458, 7159240, 
   5010068, 4317364, 2663378, 6705802, 4855975, 7946292, 676590, 7044481, 5152541, 1714295, 2453983, 1460718, 7737789, 4795319, 2815639, 2283733,
    3602218, 3182878, 2740543, 4793971, 5269599, 2101410, 3704823, 1159875, 394148, 928749, 1095468, 4874037, 2071829, 4361428, 3241972, 2156050,
    3415069, 1759347, 7562881, 4805951, 3756790, 6444618, 6663429, 4430364, 5483103, 3192354, 556856, 3870317, 2917338, 1853806, 3345963, 1858416,
    3073009, 1277625, 5744944, 3852015, 4183372, 5157610, 5258977, 8106357, 2508980, 2028118, 1937570, 4564692, 2811291, 5396636, 7270901, 4158088, 
    1528066, 482649, 1148858, 5418153, 7814814, 169688, 2462444, 5046034, 4213992, 4892034, 1987814, 5183169, 1736313, 235407, 5130263, 3258457, 
    5801164, 1787943, 5989328, 6125690, 3482206, 4197502, 7080401, 6018354, 7062739, 2461387, 3035980, 621164, 3901472, 7153756, 2925816, 3374250, 
    1356448, 5604662, 2683270, 5601629, 4912752, 2312838, 7727142, 7921254, 348812, 8052569, 1011223, 6026202, 4561790, 6458164, 6143691, 1744507,
    1753, 6444997, 5720892, 6924527, 2660408, 6600190, 8321269, 2772600, 1182243, 87208, 636927, 4415111, 4423672, 6084020, 5095502, 4663471, 8352605,
    822541, 1009365, 5926272, 6400920, 1596822, 4423473, 4620952, 6695264, 4969849, 2678278, 4611469, 4829411, 635956, 8129971, 5925040, 4234153,
    6607829, 2192938, 6653329, 2387513, 4768667, 8111961, 5199961, 3747250, 2296099, 1239911, 4541938, 3195676, 2642980, 1254190, 8368000, 2998219,
    141835, 8291116, 2513018, 7025525, 613238, 7070156, 6161950, 7921677, 6458423, 4040196, 4908348, 2039144, 6500539, 7561656, 6201452, 6757063,
    2105286, 6006015, 6346610, 586241, 7200804, 527981, 5637006, 6903432, 1994046, 2491325, 6987258, 507927, 7192532, 7655613, 6545891, 5346675, 
    8041997, 2647994, 3009748, 5767564, 4148469, 749577, 4357667, 3980599, 2569011, 6764887, 1723229, 1665318, 2028038, 1163598, 5011144, 3994671,
    8368538, 7009900, 3020393, 3363542, 214880, 545376, 7609976, 3105558, 7277073, 508145, 7826699, 860144, 3430436, 140244, 6866265, 6195333, 
    3123762, 2358373, 6187330, 5365997, 6663603, 2926054, 7987710, 8077412, 3531229, 4405932, 4606686, 1900052, 7598542, 1054478, 7648983]

### Testvector generation
print("### Test vector generation ###")
pk, sk = ML_DSA_44.keygen()
msg = b"Your message signed by ML_DSA"
sig = ML_DSA_44.sign(sk, msg)
assert ML_DSA_44.verify(pk, msg, sig)

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

def coeffFromThreeBytes(s: bytearray) -> int or None: # type: ignore
    val = int.from_bytes(s, "little")
    val &= 0x7FFFFF
    if val < polyRing.q:
        return val
    else:
        return None

def coeffFromHalfByte(b: int, eta: int) -> int or None: # type: ignore
    if eta == 2 and b < 15:
        return 2 - (b % 5)
    elif b < 9:
        assert eta == 4
        return 4 - b
    return None

def integerToBits(x: int, alpha: int):
    return int(bin(x)[2:].zfill(alpha)[::-1], 2) ### Bit reverse

def bitsToBytes(y: int, bitlen: int) -> bytearray:
    ba = bytearray()
    y = int(bin(y)[2:].zfill(bitlen * 32 * 8)[::-1], 2)  ### Bit reverse

    for i in range(polyRing.n // 8 * bitlen):
        ba += bytes([y & 0xff])
        y >>= 8

    return ba

def bytesToBits(z: bytearray) -> List(int): # type: ignore
    list_t = []
    for x in z:
        for j in range(8):
            list_t.append((x >> j) & 0x1)

    return list_t

class polyRing:
    q = 8380417
    n = 256

    def __init__(self):
        self.coeff = [0] * polyRing.n

    def __repr__(self):
        # return str(list(map(hex, self.coeff)))
        return str(self.coeff)
    
    def __getitem__(self, index):
        return self.coeff[index]
    
    def __neg__(self):
        tmp = self.__class__()
        for i in range(self.n):
            tmp.coeff[i] = -self.coeff[i] % self.q
        return tmp
    
    def __add__(self, other):
        tmp = self.__class__()
        for i in range(self.n):
            tmp.coeff[i] = (self.coeff[i] + other.coeff[i]) % self.q ### reduction
        return tmp

    def __sub__(self, other):
        tmp = self.__class__()
        for i in range(self.n):
            tmp.coeff[i] = (self.coeff[i] - other.coeff[i]) % self.q ### reduction
        return tmp

    def __lshift__(self, shift: int):
        tmp = self.__class__()
        for i in range(self.n):
            tmp.coeff[i] = self.coeff[i] << shift
        return tmp
    
    def __matmul__(self, other: polyRing) -> polyRing:
        ### Element-wise mult
        tmp = self.__class__()
        for i in range(self.n):
            tmp.coeff[i] = self.coeff[i] * other.coeff[i] % self.q

        return tmp

    def mod_pm(self) -> polyRing:
        tmp = self.__class__()
        half = (self.q-1)//2
        tmp.coeff = [elem - self.q if elem > half else elem for elem in self]
        return tmp
    
    @classmethod
    def rejNTTPoly(cls, rho: bytearray) -> Tq: # type: ignore
        G = hashlib.shake_128()
        G.update(rho)
        s = G.digest(10000)

        poly = cls()
        i = 0
        for j in range(256):
            coeff = None
            while coeff is None:
                coeff = coeffFromThreeBytes(s[i:i+3])
                i += 3
            poly.coeff[j] = coeff

        return poly

    @classmethod
    def rejBoundedPoly(cls, rho: bytearray, eta: int):
        H = hashlib.shake_256()
        H.update(rho)
        s = H.digest(256)

        poly = cls()
        i = 0
        j = 0
        while j < cls.n:
            z = s[i:i+1][0]
            i += 1
            z0 = coeffFromHalfByte(z & 0xf, eta)
            z1 = coeffFromHalfByte(z >> 4, eta)

            if z0 is not None:
                poly.coeff[j] = z0
                j += 1

            if z1 is not None and j < cls.n:
                poly.coeff[j] = z1
                j += 1

        return poly

    @classmethod
    def sampleInBall(cls, rho: bytearray, tau: int):
        H = hashlib.shake_256()
        H.update(rho)
        s = H.digest(1024)
        h = bytesToBits(s[:8])

        poly = cls()
        xof_pos = 8
        for i in range(256 - tau, 256, 1):
            j = s[xof_pos]
            xof_pos +=1
            while j > i:
                j = s[xof_pos]
                xof_pos +=1
            poly.coeff[i] = poly.coeff[j]
            poly.coeff[j] = (-1)**h[i-256+tau]

        return poly

    def power2round(self) -> (polyRing, polyRing): # type: ignore
        r0_poly = self.__class__()
        r1_poly = self.__class__()
        mask = 2**13 - 1
        half = 2**12
        """
            - n / 2  <= r <= n / 2
        """
        for i in range(self.n):
            rp = self.coeff[i] % self.q
            r0 = rp & mask
            r0 = r0 - 2**13 if r0 > half else r0
            r1 = (rp - r0) >> 13
            r0_poly.coeff[i] = r0
            r1_poly.coeff[i] = r1

        return r1_poly, r0_poly

    def decompose(self, alpha: int) -> (polyRing, polyRing): # type: ignore
        r0_poly = self.__class__()
        r1_poly = self.__class__()

        for i in range(self.n):
            rp = self.coeff[i] % self.q
            r0 = rp % (2*alpha)
            r0 = r0 - 2*alpha if r0 > alpha else r0
            if rp - r0 == self.q - 1:
                r1 = 0
                r0 = r0 - 1
            else:
                r1 = (rp - r0) // (alpha << 1)
            r0_poly.coeff[i] = r0
            r1_poly.coeff[i] = r1

        return r1_poly, r0_poly
    
    def highBits(self, alpha: int) -> polyRing:
        high, _ = self.decompose(alpha)
        return high

    def lowBits(self, alpha: int) -> polyRing:
        _, low = self.decompose(alpha)
        return low
    
    @staticmethod
    def makeHint(z: polyRing, r: polyRing, alpha: int) -> polyRing:
        poly = polyRing()
        poly.coeff = [int(r1 != v1) for r1, v1 in zip(r.highBits(alpha), (r+z).highBits(alpha))]
        return poly

    @staticmethod
    def useHint(h: polyRing, r: polyRing, alpha: int) -> polyRing:
        poly = polyRing()
        m = (polyRing.q - 1)//(2*alpha)
        r1, r0 = r.decompose(alpha)
        
        for i in range(256):
            if h[i] == 1:
                if r0[i] > 0:
                    poly.coeff[i] = (r1[i] + 1) % m
                else:
                    poly.coeff[i] = (r1[i] - 1) % m
            elif h[i] == 0:
                poly.coeff[i] = r1[i]
            else:
                raise ValueError("h is not 1/0")

        return poly
    
    def simpleBitPack(self, bitlen: int) -> bytearray:
        ### bitlen is bitlen b
        z = 0
        ### 実は逆順の係数を並べてリトルエンディアンにするだけでいい
        for i in range(polyRing.n):
            z <<= bitlen
            z |= integerToBits(self.coeff[i], bitlen)

        return bitsToBytes(z, bitlen)

    @staticmethod
    def simplBitUnpack(ba: bytearray, bitlen: int) -> polyRing:
        ### bitlen is bitlen b
        z = int.from_bytes(ba, 'little')
        poly = polyRing()
        mask = 2**bitlen - 1

        for i in range(256):
            poly.coeff[i] = z & mask
            z >>= bitlen

        return poly
    
    def bitPack(self, b: int, bitlen: int) -> bytearray:
        z = 0
        for coeff in reversed(self.coeff):
            z <<= bitlen
            z |= b - coeff

        return z.to_bytes((bitlen*256)//8, 'little')
    
    @staticmethod
    def bitUnpack(ba: bytearray, b: int, bitlen: int) -> polyRing:
        z = int.from_bytes(ba, 'little')
        poly = polyRing()
        mask = 2**bitlen - 1

        for i in range(256):
            coeff = b - (z & mask)
            poly.coeff[i] = coeff
            z >>= bitlen

        return poly
    
    def ntt(self):
        ntt = polyRing()
        ### Straightforward
        for i in range(self.n):
            br = int(f'{i:08b}'[::-1], 2)
            zeta = pow(1753, 2*br+1, self.q)
            x = 1
            sum = 0
            for j in range(self.n):
                sum += self.coeff[j] * x % self.q
                x = x * zeta
            ntt.coeff[i] = sum % self.q
        
        return ntt

    def intt(self):
        ### FFT version
        poly_out = copy.deepcopy(self)
        m = 256
        for len in [1, 2, 4, 8, 16, 32, 64, 128]:
            for start in range(0, 256, 2*len):
                m -= 1
                br = int(f'{m:08b}'[::-1], 2)
                zeta = -pow(1753, br, self.q)
                for j in range(start, start + len):
                    t = poly_out.coeff[j]
                    poly_out.coeff[j] = t + poly_out.coeff[j+len]
                    poly_out.coeff[j+len] = t - poly_out.coeff[j+len]
                    poly_out.coeff[j+len] = zeta * poly_out.coeff[j+len]

        for i in range(self.n):
            poly_out.coeff[i] *= 8347681
            poly_out.coeff[i] %= self.q
        
        return poly_out




class my_ml_dsa:
    def __init__(self):
        ### ML_DSA_44
        self.q = 8380417
        self.d = 13  # number of bits dropped from t
        self.tau = 39  # number of ±1 in c
        self.gamma_1 = 131072  # coefficient range of y: 2^17
        self.bitlen_gamma_1 = 17
        self.gamma_2 = 95232  # low order rounding range: (q-1)/88
        self.bitlen_w1 = len(bin(((self.q-1) // (2*self.gamma_2)) - 1)) - 2 # 6 or 5
        self.k = 4  # Dimensions of A = (k, l)
        self.l = 4  # Dimensions of A = (k, l)
        self.eta = 2  # Private key range
        self.omega = 80  # Max number of ones in hint
        self.beta = self.tau * self.eta
        self.c_tilde_bytes = 32

    def expandA(self, rho: bytearray) -> List(polyRing): # type: ignore
        A_hat = []
        for r in range(self.k):
            list_t =[]
            for s in range(self.l):
                rho_p = rho + s.to_bytes(1, 'little') + r.to_bytes(1, 'little')
                list_t.append(polyRing.rejNTTPoly(rho_p))
            A_hat.append(list_t)

        return A_hat

    def expandS(self, rho: bytearray) -> List(polyRing): # type: ignore
        s1 = [polyRing.rejBoundedPoly(rho + int.to_bytes(r, 2, 'little'), self.eta) for r in range(self.l)]
        s2 = [polyRing.rejBoundedPoly(rho + int.to_bytes(r + self.l, 2, 'little'), self.eta) for r in range(self.k)]

        return s1, s2

    def expandMask(self, rho: bytearray, mu: bytearray) -> List(polyRing): # type: ignore
        c = self.bitlen_gamma_1 + 1 ### 18 or 20
        y = []
        for r in range(self.l):
            rho_p = rho + int.to_bytes(mu + r, 2, 'little')
            v = hash_H(rho_p, 32*c)
            y.append(polyRing.bitUnpack(v, self.gamma_1, c))

        return y

    def w1Encode(self, w1: List(polyRing)) -> bytearray: # type: ignore
        w1_tilde = bytearray()
        for i in range(self.k):
            w1_tilde += w1[i].simpleBitPack(self.bitlen_w1)
        return w1_tilde
    
    def pkEncode(self, rho: bytearray, t1: List(polyRing)) -> bytearray: # type: ignore
        pk = rho
        for i in range(self.k):
            pk += t1[i].simpleBitPack(10)
        
        return pk
    
    def pkDecode(self, pk: bytearray):
        rho = pk[:32]
        t1 = [polyRing.simplBitUnpack(pk[32+i*320:], 10) for i in range(self.k)]

        return rho, t1

    
    def skEncode(self, rho: bytearray, K: bytearray, tr: bytearray, s1: List(polyRing), s2: List(polyRing), t0: List(List(int))) -> bytearray: # type: ignore
        sk = rho + K + tr
        bitlen = 3 if self.eta == 2 else 4 if self.eta == 4 else None

        for i in range(self.l):
            sk += s1[i].bitPack(self.eta, bitlen)
        for i in range(self.k):
            sk += s2[i].bitPack(self.eta, bitlen)
        for i in range(self.k):
            sk += t0[i].bitPack(2**12, 13)

        return sk

    def skDecode(self, sk: bytearray):
        rho, K, tr = sk[:32], sk[32:64], sk[64:128]

        s_bytes = 96 if self.eta == 2 else 128
        s1_offset = 128
        s2_offset = 128 + self.l * s_bytes
        t0_offset = 128 + 2*(self.l * s_bytes)
        bitlen = 3 if self.eta == 2 else 4 if self.eta == 4 else None

        s1 = [polyRing.bitUnpack(sk[s1_offset+i*s_bytes:], self.eta, bitlen) for i in range(self.l)]
        s2 = [polyRing.bitUnpack(sk[s2_offset+i*s_bytes:], self.eta, bitlen) for i in range(self.k)]
        t0 = [polyRing.bitUnpack(sk[t0_offset+i*416:], 2**12, 13) for i in range(self.k)]

        return rho, K, tr, s1, s2, t0

    def sigEncode(self, c_tilde: bytearray, z: List(polyRing), h: List(polyRing)) -> bytearray: # type: ignore
        sigma = c_tilde

        for i in range(self.l):
            sigma += z[i].bitPack(self.gamma_1, self.bitlen_gamma_1 + 1)
        
        ### HintBitaPack
        index = 0
        ba_index = bytearray()
        for poly in h:
            for i in range(256):
                if poly.coeff[i]:
                    sigma += i.to_bytes(1, 'big')
                    index += 1
            ba_index += index.to_bytes(1, 'big')
        sigma += bytearray([0]*(self.omega - index)) + ba_index

        return sigma
    
    def sigDecode(self, sig: bytearray):
        c_tilde = sig[:self.c_tilde_bytes]
        poly_bytes = (self.bitlen_gamma_1 + 1)*32

        z = [polyRing.bitUnpack(sig[self.c_tilde_bytes + i*poly_bytes:], self.gamma_1, self.bitlen_gamma_1 + 1) for i in range(self.l)]
        
        one_indices = sig[-(self.k+self.omega):]
        ### HinBitUnpack
        h = []
        index = 0
        for stop_index in sig[-self.k:]:
            poly = polyRing()
            while index < stop_index:
                poly.coeff[one_indices[index]] = 1
                index += 1
            h.append(poly)

        return c_tilde, z, h

    def matrix_vector_mult(self, mat: List(List(polyRing)), vec: List(polyRing)) -> List(polyRing): # type: ignore
        prod = [polyRing() for x in range(self.k)]
        for i in range(self.k):
            for j in range(self.l):
                prod[i] = prod[i] + (mat[i][j] @ vec[j])
        return prod

    def scalarVectorNTT(self, scalar: polyRing, vec: List(polyRing)) -> List(polyRing): # type: ignore
        prod = []
        for poly in vec:
            prod.append(scalar @ poly)
        return prod

    def inf_norm(self, vec: List(polyRing)) -> int: # type: ignore
        list_elem = []
        for poly in vec:
            list_elem += poly.coeff

        half = (self.q - 1)//2
        list_elem = [abs(elem - self.q) if elem > half else elem for elem in list_elem]

        return max(list_elem)

    def _keygen_internal(self, xi: int):
        hash_in = i2b(xi) + self.k.to_bytes(1, 'little') + self.l.to_bytes(1, 'little')
        seeds = hash_H(hash_in, 128)

        rho, rho_p, K = seeds[:32], seeds[32:96], seeds[96:]
        assert b2i(rho) == tv_rho, f'{rho} != {tv_rho:x}'
        assert b2i(rho_p) == tv_rho_p, f'{rho_p} != {tv_rho_p:x}' 
        assert b2i(K) == tv_K, f'{K} != {tv_K:x}' 

        A_hat = self.expandA(rho)
        s1, s2 = self.expandS(rho_p)
        s1_hat = [x.ntt() for x in s1]

        As = self.matrix_vector_mult(A_hat, s1_hat)
        t = [As[i].intt() + s2[i] for i in range(self.k)]

        ### Power2Round
        t0 = []
        t1 = []
        for poly in t:
            t1_poly, t0_poly = poly.power2round()
            t1.append(t1_poly)
            t0.append(t0_poly)

        pk = self.pkEncode(rho, t1)
        tr = hash_H(pk, 64)
        sk = self.skEncode(rho, K, tr, s1, s2, t0)
        assert b2i(pk) == tv_pk, f'{pk} != {tv_pk:x}' 
        assert b2i(tr) == tv_tr, f'{tr} != {tv_tr:x}' 
        assert b2i(sk) == tv_sk, f'{sk} != {tv_sk:x}' 

        return pk, sk

    def _sign_internal(self, sk: bytearray, Mp: bytearray, rnd: int):
        rho, K, tr, s1, s2, t0 = self.skDecode(sk)
        s1_hat = [x.ntt() for x in s1]
        s2_hat = [x.ntt() for x in s2]
        t0_hat = [x.ntt() for x in t0]
        A_hat = self.expandA(rho)
        mu = hash_H(tr + Mp, 64)
        rho_pp = hash_H(K + i2b(rnd) + mu, 64)

        kappa = 0
        z, h = None, None
        while z is None and h is None:
            y = self.expandMask(rho_pp, kappa)
            y_hat = [elem.ntt() for elem in y]
            w_hat = self.matrix_vector_mult(A_hat, y_hat)
            w = [elem.intt() for elem in w_hat]

            ### HighBits
            w1 = [poly.highBits(self.gamma_2) for poly in w]

            c_tilde = hash_H(mu + self.w1Encode(w1), self.c_tilde_bytes)
            c = polyRing.sampleInBall(c_tilde, self.tau)
            c_hat = c.ntt()

            cs1 = [x.intt() for x in self.scalarVectorNTT(c_hat, s1_hat)]
            cs2 = [x.intt() for x in self.scalarVectorNTT(c_hat, s2_hat)]
            z = [y[i] + cs1[i] for i in range(self.l)]

            wcs = [w[i] - cs2[i] for i in range(self.k)]
            ### LowBits
            r0 = [poly.lowBits(self.gamma_2) for poly in wcs]


            ### Validity check
            if self.inf_norm(z) >= self.gamma_1 - self.beta or self.inf_norm(r0) >= self.gamma_2 - self.beta:
                z, h = None, None
            else:
                ct0 = [x.intt() for x in self.scalarVectorNTT(c_hat, t0_hat)]
                h = [polyRing.makeHint(-ct0[i], w[i] - cs2[i] + ct0[i], self.gamma_2) for i in range(self.k)]
                if self.inf_norm(ct0) > self.gamma_2 or sum([sum(poly.coeff) for poly in h]) > self.omega:
                    z, h = None, None
            kappa += self.l

        z = [poly.mod_pm() for poly in z]
        sigma = self.sigEncode(c_tilde, z, h)
        assert b2i(sigma) == tv_sig, f'{sigma} != {tv_sig:x}' 

        return sigma
    
    def _verify_internal(self, pk: bytearray, Mp: bytearray, sig: bytearray) -> bool:
        rho, t1 = self.pkDecode(pk)
        c_tilde, z, h = self.sigDecode(sig)

        A_hat = self.expandA(rho)
        tr = hash_H(pk, 64)
        mu = hash_H(tr + Mp, 64)
        c = polyRing.sampleInBall(c_tilde, self.tau)

        z_hat = [elem.ntt() for elem in z]
        c_hat = c.ntt()
        t1_hat = [(elem << self.d).ntt() for elem in t1]

        Az_hat = self.matrix_vector_mult(A_hat, z_hat)
        ct1_hat = self.scalarVectorNTT(c_hat, t1_hat)
        wp_approx = [(Az_hat[i] - ct1_hat[i]).intt() for i in range(self.k)]

        w1_p = [polyRing.useHint(h[i], wp_approx[i], self.gamma_2) for i in range(self.k)]
        c_tilde_p = hash_H(mu + self.w1Encode(w1_p), self.c_tilde_bytes)

        return self.inf_norm(z) < self.gamma_1 - self.beta and c_tilde == c_tilde_p


    def keygen(self):
        xi = tv_xi
        return self._keygen_internal(xi)

    def sign(self, sk: bytearray, M: bytearray, ctx: bytearray = b"") -> bytearray:
        if len(ctx) > 255:
            raise ValueError("ctx length > 255")
        
        rnd = tv_rnd

        Mp = bytes([0, len(ctx)]) + ctx + M
        sig = self._sign_internal(sk, Mp, rnd)

        return sig

    def verify(self, pk: bytearray, M: bytearray, sig: bytearray, ctx: bytearray = b"") -> bytearray:
        if len(ctx) > 255:
            raise ValueError("ctx length > 255")
        
        Mp = bytes([0, len(ctx)]) + ctx + M
        return self._verify_internal(pk, Mp, sig)

inst = my_ml_dsa()
pk, sk = inst.keygen()
sig = inst.sign(sk, b"Your message signed by ML_DSA")
res = inst.verify(pk, b"Your message signed by ML_DSA", sig)
print(res)

