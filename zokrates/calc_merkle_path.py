import subprocess
import hashlib
import binascii

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

'''
Generates merkle tree path, given a
''' 
def calc_merkle(index, txids):
    # TODO


root_LE = flip32Bytes(root)
left_LE = flip32Bytes(left)
right_LE = flip32Bytes(right)
print(root_LE)
print(double_sha256(left_LE+right_LE))