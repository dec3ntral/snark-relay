

import subprocess
import hashlib
import binascii

def halve_sha256(string):
    firstpart, secondpart = string[::2], string[1::2] 
    return (firstpart, secondpart)

def pack_ints(nhigh, nlow, bits):
    return (nhigh << bits) + nlow

def unpack_int(n, bits): 
    mask = (2 ** bits - 1) 
    n_low = n & mask 
    n_high = (n >> bits) & mask 
    return [n_high, n_low]

def str_to_int(string):
    return int(string,16)

def construct_arguments(root, proof, txid):
    arguments = []
    # root
    arguments.extend(
        unpack_int(str_to_int(root), 128)
    )
    arguments.extend(
                unpack_int(str_to_int(proof[0][0]), 128)
            )
    arguments.append(proof[0][1])

    # proof - max length = 33
    '''
    for i in range(0,11):
        if i < len(proof):
            arguments.extend(
                unpack_int(str_to_int(proof[i][0]), 128)
            )
            arguments.append(proof[i][1])
        else:
            arguments.extend([0, 0, 0])

    '''
    # txid to be proven
    arguments.extend(
        unpack_int(str_to_int(txid), 128)
    )
    return [str(v) for v in arguments]



# root 
root = "cb0c453a3f456ed6eeabb1981cde22fda18233fd2665011705b21c5fbe154f17"
#left: 
left = "1ff36423b25dd0d0f3604190476d64ac593891d69d9595dceabf6e400e055deb"
#right:
right = "b24d0b39bb06e8405d3658e9b74a6efb2c7e8898fa2205a30a19a390f12d816b"


def dblShaFlip(header):
    first_hash = hashlib.sha256(binascii.unhexlify(header)).hexdigest()
    second_hash = hashlib.sha256(binascii.unhexlify(first_hash)).hexdigest()
    return flip32Bytes(second_hash)

def flip32Bytes(b32):
    byteSize = 2
    chunks = [ b32[i:i+byteSize] for i in range(0, len(b32), byteSize) ]
    reversed_chunks = chunks[::-1]
    return ('').join(reversed_chunks)

def double_sha256(data):
    hash = hashlib.sha256(binascii.unhexlify(data)).hexdigest()
    hash2 = hashlib.sha256(binascii.unhexlify(hash)).hexdigest()
    return hash2

root_LE = flip32Bytes(root)
left_LE = flip32Bytes(left)
right_LE = flip32Bytes(right)
print(root_LE)
print(double_sha256(left_LE+right_LE))

# Zokrates function signature:
# (field[2] root, field[33] proof, field[2] tx_hash)

# proof: [(left-hash, 0), (right-hash, 1)]

arguments = construct_arguments(root_LE, [(left_LE, 1)], right_LE)
print(arguments)
subprocess.call(["zokrates", "compute-witness", "-a", *arguments])

