from bitcoinrpc.authproxy import AuthServiceProxy, JSONRPCException
import pprint
import merkletools
import base58
import binascii

# 000000000000007ecdf2fc6af2f7a566161c35c4a185fcb308ab88135aedf047

pp = pprint.PrettyPrinter(indent=4)

rpc_user="user"
rpc_password="password"




rpc_connection = AuthServiceProxy("http://%s:%s@127.0.0.1:8332"%(rpc_user, rpc_password))
block =  rpc_connection.getblock("000000000000007ecdf2fc6af2f7a566161c35c4a185fcb308ab88135aedf047")

merkleProof = rpc_connection.gettxoutproof(["938a1e3dfc13d96ca29e64d15a059f06830d8289b6c16952709a4c6ab01098de"],"000000000000007ecdf2fc6af2f7a566161c35c4a185fcb308ab88135aedf047")

# print(len(merkleProof))
print(binascii.unhexlify(merkleProof))

