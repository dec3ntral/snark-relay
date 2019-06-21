// This file is LGPL3 Licensed

/**
 * @title Elliptic curve operations on twist points for alt_bn128
 * @author Mustafa Al-Bassam (mus@musalbas.com)
 */
library BN256G2 {
    uint256 internal constant FIELD_MODULUS = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47;
    uint256 internal constant TWISTBX = 0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5;
    uint256 internal constant TWISTBY = 0x9713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2;
    uint internal constant PTXX = 0;
    uint internal constant PTXY = 1;
    uint internal constant PTYX = 2;
    uint internal constant PTYY = 3;
    uint internal constant PTZX = 4;
    uint internal constant PTZY = 5;

    /**
     * @notice Add two twist points
     * @param pt1xx Coefficient 1 of x on point 1
     * @param pt1xy Coefficient 2 of x on point 1
     * @param pt1yx Coefficient 1 of y on point 1
     * @param pt1yy Coefficient 2 of y on point 1
     * @param pt2xx Coefficient 1 of x on point 2
     * @param pt2xy Coefficient 2 of x on point 2
     * @param pt2yx Coefficient 1 of y on point 2
     * @param pt2yy Coefficient 2 of y on point 2
     * @return (pt3xx, pt3xy, pt3yx, pt3yy)
     */
    function ECTwistAdd(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy
    ) public pure returns (
        uint256, uint256,
        uint256, uint256
    ) {
        if (
            pt1xx == 0 && pt1xy == 0 &&
            pt1yx == 0 && pt1yy == 0
        ) {
            if (!(
                pt2xx == 0 && pt2xy == 0 &&
                pt2yx == 0 && pt2yy == 0
            )) {
                assert(_isOnCurve(
                    pt2xx, pt2xy,
                    pt2yx, pt2yy
                ));
            }
            return (
                pt2xx, pt2xy,
                pt2yx, pt2yy
            );
        } else if (
            pt2xx == 0 && pt2xy == 0 &&
            pt2yx == 0 && pt2yy == 0
        ) {
            assert(_isOnCurve(
                pt1xx, pt1xy,
                pt1yx, pt1yy
            ));
            return (
                pt1xx, pt1xy,
                pt1yx, pt1yy
            );
        }

        assert(_isOnCurve(
            pt1xx, pt1xy,
            pt1yx, pt1yy
        ));
        assert(_isOnCurve(
            pt2xx, pt2xy,
            pt2yx, pt2yy
        ));

        uint256[6] memory pt3 = _ECTwistAddJacobian(
            pt1xx, pt1xy,
            pt1yx, pt1yy,
            1,     0,
            pt2xx, pt2xy,
            pt2yx, pt2yy,
            1,     0
        );

        return _fromJacobian(
            pt3[PTXX], pt3[PTXY],
            pt3[PTYX], pt3[PTYY],
            pt3[PTZX], pt3[PTZY]
        );
    }

    /**
     * @notice Multiply a twist point by a scalar
     * @param s     Scalar to multiply by
     * @param pt1xx Coefficient 1 of x
     * @param pt1xy Coefficient 2 of x
     * @param pt1yx Coefficient 1 of y
     * @param pt1yy Coefficient 2 of y
     * @return (pt2xx, pt2xy, pt2yx, pt2yy)
     */
    function ECTwistMul(
        uint256 s,
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy
    ) public pure returns (
        uint256, uint256,
        uint256, uint256
    ) {
        uint256 pt1zx = 1;
        if (
            pt1xx == 0 && pt1xy == 0 &&
            pt1yx == 0 && pt1yy == 0
        ) {
            pt1xx = 1;
            pt1yx = 1;
            pt1zx = 0;
        } else {
            assert(_isOnCurve(
                pt1xx, pt1xy,
                pt1yx, pt1yy
            ));
        }

        uint256[6] memory pt2 = _ECTwistMulJacobian(
            s,
            pt1xx, pt1xy,
            pt1yx, pt1yy,
            pt1zx, 0
        );

        return _fromJacobian(
            pt2[PTXX], pt2[PTXY],
            pt2[PTYX], pt2[PTYY],
            pt2[PTZX], pt2[PTZY]
        );
    }

    /**
     * @notice Get the field modulus
     * @return The field modulus
     */
    function GetFieldModulus() public pure returns (uint256) {
        return FIELD_MODULUS;
    }

    function submod(uint256 a, uint256 b, uint256 n) internal pure returns (uint256) {
        return addmod(a, n - b, n);
    }

    function _FQ2Mul(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns(uint256, uint256) {
        return (
            submod(mulmod(xx, yx, FIELD_MODULUS), mulmod(xy, yy, FIELD_MODULUS), FIELD_MODULUS),
            addmod(mulmod(xx, yy, FIELD_MODULUS), mulmod(xy, yx, FIELD_MODULUS), FIELD_MODULUS)
        );
    }

    function _FQ2Muc(
        uint256 xx, uint256 xy,
        uint256 c
    ) internal pure returns(uint256, uint256) {
        return (
            mulmod(xx, c, FIELD_MODULUS),
            mulmod(xy, c, FIELD_MODULUS)
        );
    }

    function _FQ2Add(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns(uint256, uint256) {
        return (
            addmod(xx, yx, FIELD_MODULUS),
            addmod(xy, yy, FIELD_MODULUS)
        );
    }

    function _FQ2Sub(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns(uint256 rx, uint256 ry) {
        return (
            submod(xx, yx, FIELD_MODULUS),
            submod(xy, yy, FIELD_MODULUS)
        );
    }

    function _FQ2Div(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns(uint256, uint256) {
        (yx, yy) = _FQ2Inv(yx, yy);
        return _FQ2Mul(xx, xy, yx, yy);
    }

    function _FQ2Inv(uint256 x, uint256 y) internal pure returns(uint256, uint256) {
        uint256 inv = _modInv(addmod(mulmod(y, y, FIELD_MODULUS), mulmod(x, x, FIELD_MODULUS), FIELD_MODULUS), FIELD_MODULUS);
        return (
            mulmod(x, inv, FIELD_MODULUS),
            FIELD_MODULUS - mulmod(y, inv, FIELD_MODULUS)
        );
    }

    function _isOnCurve(
        uint256 xx, uint256 xy,
        uint256 yx, uint256 yy
    ) internal pure returns (bool) {
        uint256 yyx;
        uint256 yyy;
        uint256 xxxx;
        uint256 xxxy;
        (yyx, yyy) = _FQ2Mul(yx, yy, yx, yy);
        (xxxx, xxxy) = _FQ2Mul(xx, xy, xx, xy);
        (xxxx, xxxy) = _FQ2Mul(xxxx, xxxy, xx, xy);
        (yyx, yyy) = _FQ2Sub(yyx, yyy, xxxx, xxxy);
        (yyx, yyy) = _FQ2Sub(yyx, yyy, TWISTBX, TWISTBY);
        return yyx == 0 && yyy == 0;
    }

    function _modInv(uint256 a, uint256 n) internal pure returns(uint256 t) {
        t = 0;
        uint256 newT = 1;
        uint256 r = n;
        uint256 newR = a;
        uint256 q;
        while (newR != 0) {
            q = r / newR;
            (t, newT) = (newT, submod(t, mulmod(q, newT, n), n));
            (r, newR) = (newR, r - q * newR);
        }
    }

    function _fromJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal pure returns (
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy
    ) {
        uint256 invzx;
        uint256 invzy;
        (invzx, invzy) = _FQ2Inv(pt1zx, pt1zy);
        (pt2xx, pt2xy) = _FQ2Mul(pt1xx, pt1xy, invzx, invzy);
        (pt2yx, pt2yy) = _FQ2Mul(pt1yx, pt1yy, invzx, invzy);
    }

    function _ECTwistAddJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy,
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy,
        uint256 pt2zx, uint256 pt2zy) internal pure returns (uint256[6] memory pt3) {
            if (pt1zx == 0 && pt1zy == 0) {
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    pt2xx, pt2xy,
                    pt2yx, pt2yy,
                    pt2zx, pt2zy
                );
                return pt3;
            } else if (pt2zx == 0 && pt2zy == 0) {
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    pt1xx, pt1xy,
                    pt1yx, pt1yy,
                    pt1zx, pt1zy
                );
                return pt3;
            }

            (pt2yx,     pt2yy)     = _FQ2Mul(pt2yx, pt2yy, pt1zx, pt1zy); // U1 = y2 * z1
            (pt3[PTYX], pt3[PTYY]) = _FQ2Mul(pt1yx, pt1yy, pt2zx, pt2zy); // U2 = y1 * z2
            (pt2xx,     pt2xy)     = _FQ2Mul(pt2xx, pt2xy, pt1zx, pt1zy); // V1 = x2 * z1
            (pt3[PTZX], pt3[PTZY]) = _FQ2Mul(pt1xx, pt1xy, pt2zx, pt2zy); // V2 = x1 * z2

            if (pt2xx == pt3[PTZX] && pt2xy == pt3[PTZY]) {
                if (pt2yx == pt3[PTYX] && pt2yy == pt3[PTYY]) {
                    (
                        pt3[PTXX], pt3[PTXY],
                        pt3[PTYX], pt3[PTYY],
                        pt3[PTZX], pt3[PTZY]
                    ) = _ECTwistDoubleJacobian(pt1xx, pt1xy, pt1yx, pt1yy, pt1zx, pt1zy);
                    return pt3;
                }
                (
                    pt3[PTXX], pt3[PTXY],
                    pt3[PTYX], pt3[PTYY],
                    pt3[PTZX], pt3[PTZY]
                ) = (
                    1, 0,
                    1, 0,
                    0, 0
                );
                return pt3;
            }

            (pt2zx,     pt2zy)     = _FQ2Mul(pt1zx, pt1zy, pt2zx,     pt2zy);     // W = z1 * z2
            (pt1xx,     pt1xy)     = _FQ2Sub(pt2yx, pt2yy, pt3[PTYX], pt3[PTYY]); // U = U1 - U2
            (pt1yx,     pt1yy)     = _FQ2Sub(pt2xx, pt2xy, pt3[PTZX], pt3[PTZY]); // V = V1 - V2
            (pt1zx,     pt1zy)     = _FQ2Mul(pt1yx, pt1yy, pt1yx,     pt1yy);     // V_squared = V * V
            (pt2yx,     pt2yy)     = _FQ2Mul(pt1zx, pt1zy, pt3[PTZX], pt3[PTZY]); // V_squared_times_V2 = V_squared * V2
            (pt1zx,     pt1zy)     = _FQ2Mul(pt1zx, pt1zy, pt1yx,     pt1yy);     // V_cubed = V * V_squared
            (pt3[PTZX], pt3[PTZY]) = _FQ2Mul(pt1zx, pt1zy, pt2zx,     pt2zy);     // newz = V_cubed * W
            (pt2xx,     pt2xy)     = _FQ2Mul(pt1xx, pt1xy, pt1xx,     pt1xy);     // U * U
            (pt2xx,     pt2xy)     = _FQ2Mul(pt2xx, pt2xy, pt2zx,     pt2zy);     // U * U * W
            (pt2xx,     pt2xy)     = _FQ2Sub(pt2xx, pt2xy, pt1zx,     pt1zy);     // U * U * W - V_cubed
            (pt2zx,     pt2zy)     = _FQ2Muc(pt2yx, pt2yy, 2);                    // 2 * V_squared_times_V2
            (pt2xx,     pt2xy)     = _FQ2Sub(pt2xx, pt2xy, pt2zx,     pt2zy);     // A = U * U * W - V_cubed - 2 * V_squared_times_V2
            (pt3[PTXX], pt3[PTXY]) = _FQ2Mul(pt1yx, pt1yy, pt2xx,     pt2xy);     // newx = V * A
            (pt1yx,     pt1yy)     = _FQ2Sub(pt2yx, pt2yy, pt2xx,     pt2xy);     // V_squared_times_V2 - A
            (pt1yx,     pt1yy)     = _FQ2Mul(pt1xx, pt1xy, pt1yx,     pt1yy);     // U * (V_squared_times_V2 - A)
            (pt1xx,     pt1xy)     = _FQ2Mul(pt1zx, pt1zy, pt3[PTYX], pt3[PTYY]); // V_cubed * U2
            (pt3[PTYX], pt3[PTYY]) = _FQ2Sub(pt1yx, pt1yy, pt1xx,     pt1xy);     // newy = U * (V_squared_times_V2 - A) - V_cubed * U2
    }

    function _ECTwistDoubleJacobian(
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal pure returns(
        uint256 pt2xx, uint256 pt2xy,
        uint256 pt2yx, uint256 pt2yy,
        uint256 pt2zx, uint256 pt2zy
    ) {
        (pt2xx, pt2xy) = _FQ2Muc(pt1xx, pt1xy, 3);            // 3 * x
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1xx, pt1xy); // W = 3 * x * x
        (pt1zx, pt1zy) = _FQ2Mul(pt1yx, pt1yy, pt1zx, pt1zy); // S = y * z
        (pt2yx, pt2yy) = _FQ2Mul(pt1xx, pt1xy, pt1yx, pt1yy); // x * y
        (pt2yx, pt2yy) = _FQ2Mul(pt2yx, pt2yy, pt1zx, pt1zy); // B = x * y * S
        (pt1xx, pt1xy) = _FQ2Mul(pt2xx, pt2xy, pt2xx, pt2xy); // W * W
        (pt2zx, pt2zy) = _FQ2Muc(pt2yx, pt2yy, 8);            // 8 * B
        (pt1xx, pt1xy) = _FQ2Sub(pt1xx, pt1xy, pt2zx, pt2zy); // H = W * W - 8 * B
        (pt2zx, pt2zy) = _FQ2Mul(pt1zx, pt1zy, pt1zx, pt1zy); // S_squared = S * S
        (pt2yx, pt2yy) = _FQ2Muc(pt2yx, pt2yy, 4);            // 4 * B
        (pt2yx, pt2yy) = _FQ2Sub(pt2yx, pt2yy, pt1xx, pt1xy); // 4 * B - H
        (pt2yx, pt2yy) = _FQ2Mul(pt2yx, pt2yy, pt2xx, pt2xy); // W * (4 * B - H)
        (pt2xx, pt2xy) = _FQ2Muc(pt1yx, pt1yy, 8);            // 8 * y
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1yx, pt1yy); // 8 * y * y
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt2zx, pt2zy); // 8 * y * y * S_squared
        (pt2yx, pt2yy) = _FQ2Sub(pt2yx, pt2yy, pt2xx, pt2xy); // newy = W * (4 * B - H) - 8 * y * y * S_squared
        (pt2xx, pt2xy) = _FQ2Muc(pt1xx, pt1xy, 2);            // 2 * H
        (pt2xx, pt2xy) = _FQ2Mul(pt2xx, pt2xy, pt1zx, pt1zy); // newx = 2 * H * S
        (pt2zx, pt2zy) = _FQ2Mul(pt1zx, pt1zy, pt2zx, pt2zy); // S * S_squared
        (pt2zx, pt2zy) = _FQ2Muc(pt2zx, pt2zy, 8);            // newz = 8 * S * S_squared
    }

    function _ECTwistMulJacobian(
        uint256 d,
        uint256 pt1xx, uint256 pt1xy,
        uint256 pt1yx, uint256 pt1yy,
        uint256 pt1zx, uint256 pt1zy
    ) internal pure returns(uint256[6] memory pt2) {
        while (d != 0) {
            if ((d & 1) != 0) {
                pt2 = _ECTwistAddJacobian(
                    pt2[PTXX], pt2[PTXY],
                    pt2[PTYX], pt2[PTYY],
                    pt2[PTZX], pt2[PTZY],
                    pt1xx, pt1xy,
                    pt1yx, pt1yy,
                    pt1zx, pt1zy);
            }
            (
                pt1xx, pt1xy,
                pt1yx, pt1yy,
                pt1zx, pt1zy
            ) = _ECTwistDoubleJacobian(
                pt1xx, pt1xy,
                pt1yx, pt1yy,
                pt1zx, pt1zy
            );

            d = d / 2;
        }
    }
}


// This file is MIT Licensed.
//
// Copyright 2017 Christian Reitwiessner
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

pragma solidity ^0.5.0;
library Pairing {
    struct G1Point {
        uint X;
        uint Y;
    }
    // Encoding of field elements is: X[0] * z + X[1]
    struct G2Point {
        uint[2] X;
        uint[2] Y;
    }
    /// @return the generator of G1
    function P1() pure internal returns (G1Point memory) {
        return G1Point(1, 2);
    }
    /// @return the generator of G2
    function P2() pure internal returns (G2Point memory) {
        return G2Point(
            [11559732032986387107991004021392285783925812861821192530917403151452391805634,
             10857046999023057135944570762232829481370756359578518086990519993285655852781],
            [4082367875863433681332203403145435568316851327593401208105741076214120093531,
             8495653923123431417604973247489272438418190587263600148770280649306958101930]
        );
    }
    /// @return the negation of p, i.e. p.addition(p.negate()) should be zero.
    function negate(G1Point memory p) pure internal returns (G1Point memory) {
        // The prime q in the base field F_q for G1
        uint q = 21888242871839275222246405745257275088696311157297823662689037894645226208583;
        if (p.X == 0 && p.Y == 0)
            return G1Point(0, 0);
        return G1Point(p.X, q - (p.Y % q));
    }
    /// @return the sum of two points of G1
    function addition(G1Point memory p1, G1Point memory p2) internal returns (G1Point memory r) {
        uint[4] memory input;
        input[0] = p1.X;
        input[1] = p1.Y;
        input[2] = p2.X;
        input[3] = p2.Y;
        bool success;
        assembly {
            success := call(sub(gas, 2000), 6, 0, input, 0xc0, r, 0x60)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require(success);
    }
    /// @return the sum of two points of G2
    function addition(G2Point memory p1, G2Point memory p2) internal pure returns (G2Point memory r) {
        (r.X[1], r.X[0], r.Y[1], r.Y[0]) = BN256G2.ECTwistAdd(p1.X[1],p1.X[0],p1.Y[1],p1.Y[0],p2.X[1],p2.X[0],p2.Y[1],p2.Y[0]);
    }
    /// @return the product of a point on G1 and a scalar, i.e.
    /// p == p.scalar_mul(1) and p.addition(p) == p.scalar_mul(2) for all points p.
    function scalar_mul(G1Point memory p, uint s) internal returns (G1Point memory r) {
        uint[3] memory input;
        input[0] = p.X;
        input[1] = p.Y;
        input[2] = s;
        bool success;
        assembly {
            success := call(sub(gas, 2000), 7, 0, input, 0x80, r, 0x60)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require (success);
    }
    /// @return the result of computing the pairing check
    /// e(p1[0], p2[0]) *  .... * e(p1[n], p2[n]) == 1
    /// For example pairing([P1(), P1().negate()], [P2(), P2()]) should
    /// return true.
    function pairing(G1Point[] memory p1, G2Point[] memory p2) internal returns (bool) {
        require(p1.length == p2.length);
        uint elements = p1.length;
        uint inputSize = elements * 6;
        uint[] memory input = new uint[](inputSize);
        for (uint i = 0; i < elements; i++)
        {
            input[i * 6 + 0] = p1[i].X;
            input[i * 6 + 1] = p1[i].Y;
            input[i * 6 + 2] = p2[i].X[0];
            input[i * 6 + 3] = p2[i].X[1];
            input[i * 6 + 4] = p2[i].Y[0];
            input[i * 6 + 5] = p2[i].Y[1];
        }
        uint[1] memory out;
        bool success;
        assembly {
            success := call(sub(gas, 2000), 8, 0, add(input, 0x20), mul(inputSize, 0x20), out, 0x20)
            // Use "invalid" to make gas estimation work
            switch success case 0 { invalid() }
        }
        require(success);
        return out[0] != 0;
    }
    /// Convenience method for a pairing check for two pairs.
    function pairingProd2(G1Point memory a1, G2Point memory a2, G1Point memory b1, G2Point memory b2) internal returns (bool) {
        G1Point[] memory p1 = new G1Point[](2);
        G2Point[] memory p2 = new G2Point[](2);
        p1[0] = a1;
        p1[1] = b1;
        p2[0] = a2;
        p2[1] = b2;
        return pairing(p1, p2);
    }
    /// Convenience method for a pairing check for three pairs.
    function pairingProd3(
            G1Point memory a1, G2Point memory a2,
            G1Point memory b1, G2Point memory b2,
            G1Point memory c1, G2Point memory c2
    ) internal returns (bool) {
        G1Point[] memory p1 = new G1Point[](3);
        G2Point[] memory p2 = new G2Point[](3);
        p1[0] = a1;
        p1[1] = b1;
        p1[2] = c1;
        p2[0] = a2;
        p2[1] = b2;
        p2[2] = c2;
        return pairing(p1, p2);
    }
    /// Convenience method for a pairing check for four pairs.
    function pairingProd4(
            G1Point memory a1, G2Point memory a2,
            G1Point memory b1, G2Point memory b2,
            G1Point memory c1, G2Point memory c2,
            G1Point memory d1, G2Point memory d2
    ) internal returns (bool) {
        G1Point[] memory p1 = new G1Point[](4);
        G2Point[] memory p2 = new G2Point[](4);
        p1[0] = a1;
        p1[1] = b1;
        p1[2] = c1;
        p1[3] = d1;
        p2[0] = a2;
        p2[1] = b2;
        p2[2] = c2;
        p2[3] = d2;
        return pairing(p1, p2);
    }
}

contract Verifier {
    using Pairing for *;
    struct VerifyingKey {
        Pairing.G1Point a;
        Pairing.G2Point b;
        Pairing.G2Point gamma;
        Pairing.G2Point delta;
        Pairing.G1Point[] gammaABC;
    }
    struct Proof {
        Pairing.G1Point A;
        Pairing.G2Point B;
        Pairing.G1Point C;
    }
    function verifyingKey() pure internal returns (VerifyingKey memory vk) {
        vk.a = Pairing.G1Point(uint256(0x2d73c58346cd3ae69cdd821df2f37b752d3f45accb3a4a9405468af7369680fb), uint256(0x1f018f97d4148a9190b5a710bc8c92cd8eebf49ca128d29cdec6bf4c439c7418));
        vk.b = Pairing.G2Point([uint256(0x0fc0bf5a7a70a0b2d00e767ae12c6019c148c917e85715554ecfe755fda3b01a), uint256(0x13ee3ba57781400bc7953a5f6b2d3ab52a3bb83f914ef2ed19c63834a3037219)], [uint256(0x0c11ce41932ab682ff4a7274d5521fcafd8eff4a39a4fdaa83ed5fa3fa933efb), uint256(0x28be3667568c25e029aed7c0ae1b0a05f9371b6dac9631cec867f6668fbe798a)]);
        vk.gamma = Pairing.G2Point([uint256(0x08d6e56a1e67c1c9f59bc75d6f98654a2fdb0c46fe7966aa5a2936c0492ca4e1), uint256(0x22c0183680bf17c038f16281269205b61361157b5a267fcc3d3add6b4654cfc7)], [uint256(0x1159411b9624269aa6d681016a3aa53c3477f93b4252f43f3ae14fb323656255), uint256(0x1ee64f7c6a7b11505be3f64dc045fbf8301ce81f7b8b814f063dd7a82b15588f)]);
        vk.delta = Pairing.G2Point([uint256(0x227d93acf1ca9d784af09053749da84b9c3270cf18b0626418c2a021133c95be), uint256(0x27b62bd73f07dd18742a9843f28ac7002d8cf79974849909f698bf1c02713ee0)], [uint256(0x2d7c42302e36b0c3ae81f2f724130676c9832646e3eb2a08a9f443b831b5de51), uint256(0x24081c75d5c91abcd0549f09713d91d3185e87b64e01592b8f93d70034c68bcf)]);
        vk.gammaABC = new Pairing.G1Point[](39);
        vk.gammaABC[0] = Pairing.G1Point(uint256(0x2470f86b6cff2b678e9637803e5bf3f50647befa812d40012bcbdd25778703a0), uint256(0x25e065e13c895c4ebe53c16f3513c482b9cae154bd6014080185a73a23d1f6f3));
        vk.gammaABC[1] = Pairing.G1Point(uint256(0x1211e32f29b8860775fe94f4b3c72542efa5a3ac866a99cfc2fabe403d879cb1), uint256(0x298d724ce136eb9c0ab82ca37743e87a92dcfceb22efeba96aef5672fdbc8b7c));
        vk.gammaABC[2] = Pairing.G1Point(uint256(0x094257d31c1e9f66656af7c2e013c22586104a730bafba47f29dc5480e2be2be), uint256(0x0a3aeffae5465241630e2c5cd8bf2ef46db84ecf3e8f9e213e5d340061d9d224));
        vk.gammaABC[3] = Pairing.G1Point(uint256(0x0362e3f127b3324ca73a1aed2d31e56463c23d4be762355995685a65194c3f1e), uint256(0x124f58c2395620b07c8854bdd0ffa3ade185c902f137c5ef291ebd12fa574871));
        vk.gammaABC[4] = Pairing.G1Point(uint256(0x1e71c06be506a9bceef971e977133bef7d4df8cd449224f2f034bb29727b9536), uint256(0x02e10a15499f1888067af43ccd36440ac40905f3935dd4d58342a4dbc68be982));
        vk.gammaABC[5] = Pairing.G1Point(uint256(0x0248ae5f303f33c1941745f99acc6207cf2b21d80e31437af499eca1892c4c9c), uint256(0x0b30509851eb2c63bd9f8f208fca7ce739294f61e8e5d05385b1e9a2f2c9eb9c));
        vk.gammaABC[6] = Pairing.G1Point(uint256(0x1ce6f854cf01fb1d202a1400e64678b5564a4ed4ca460ff761dd6859d37df9cc), uint256(0x0c79274db48c7c8c7bccbceb37e6110dd0bad788b1ae7c5169be21e4d6c1c341));
        vk.gammaABC[7] = Pairing.G1Point(uint256(0x26e858521d54bbc81974745202ac0941677ce38d8fe2b9789d3a22c9e87e7be6), uint256(0x22aaae7dcfffb8f6fa5a5fa77565a4c3327478b4e5a0b24e35e592c4dd2494cb));
        vk.gammaABC[8] = Pairing.G1Point(uint256(0x0a59ec716d2e164ed65f4596669ea3f109a4c30586a543ffbc66b4d422189e02), uint256(0x14ba53a43175f93029dcb6e4c385f8ee6cd97033ac4d57f3b62a54fa70bc7654));
        vk.gammaABC[9] = Pairing.G1Point(uint256(0x0a55d7f0c3dc7cf5980e7c354b7386cfd96f0dcdf7b9a7ab605877afd91bbbe5), uint256(0x18f3a03cd7cc82d11980c349cbd3cf2b69837ce5cc7da47979e593df197f0fc6));
        vk.gammaABC[10] = Pairing.G1Point(uint256(0x08439093f1cf3c13ef627142415f440f56775bf8fd5ed21baf74a8efa4cad053), uint256(0x2f9929fbe9e927698437325e3f5d11ddcc4cd222817d9a4d8305c870dab2dcb9));
        vk.gammaABC[11] = Pairing.G1Point(uint256(0x29252e1eeac83fa0a5ea2d0d1a96f97a689299a2f200b540c86efaf8c9005e21), uint256(0x14690d46f8b7a6288d06033d67e94b9ae1633e3410cd40c6381c49cb65c221e4));
        vk.gammaABC[12] = Pairing.G1Point(uint256(0x0ce2f817747c3475d4fcc48f2140d33c42a21ac48ca83486f99146d9cc602338), uint256(0x1a2a80b782ffcadecb494b26bc0b7724cfe70fa3b50c509a7f854c2f331a6887));
        vk.gammaABC[13] = Pairing.G1Point(uint256(0x138e120ea4cba5a56dbbdb318e5c3d9960516c25e0744efdf6f1cbfbdebe0d50), uint256(0x06081bde5ecb251f3911caad59eb0bf89e6f038dd4dc0692816d71de48251f1b));
        vk.gammaABC[14] = Pairing.G1Point(uint256(0x027ae4bb7c381ff5db18999fa75e2339dc8ce261142ec6fffe81dfb1758ce384), uint256(0x0790f52ba743d67ee3a1f88369af5d70d7729950e3944d621cd77e25eee9443a));
        vk.gammaABC[15] = Pairing.G1Point(uint256(0x008df65ffa2c7ef6d6ca86154d0ae39b0a2bf1f6b363b0b5eb9947f903fd2b28), uint256(0x16a27fec3c7dc30dc023a2222a7ce1526be5d0b4633478334a3f0df8ddbfb5ca));
        vk.gammaABC[16] = Pairing.G1Point(uint256(0x27407148723f7ad405131a0927b666a7e8520b5253255fa42c00d5702a2bf094), uint256(0x10b6fe183dca032afe5e8f0fd153a16558103dfa2092d109a4628db79d8b4d66));
        vk.gammaABC[17] = Pairing.G1Point(uint256(0x0ef4e76bffb9e0cda5ee872953b4530d331bbc5128d0723519fa15b6d3b20422), uint256(0x156b7da1b71eb0df73d576d2d7fcaa62316b841c8be2c4d6adfd8f4d8fbbb925));
        vk.gammaABC[18] = Pairing.G1Point(uint256(0x18f8973e393b6f3185befe3886b34cf5877081ec64e78ecd34dbc10b6b29bdf3), uint256(0x0da002359ac46f56a5238f23a1be2bb3bfbe7ffe895a2db933877a9a9b86da37));
        vk.gammaABC[19] = Pairing.G1Point(uint256(0x17db572a07fd72571a514e9d5bd9a1cdc4dc9a2b7f7642036699b23f4fdd7142), uint256(0x215fdfa481d6f09ed3869741ffa4fe9974b1d1ea48c4d55a67b8d85648569e12));
        vk.gammaABC[20] = Pairing.G1Point(uint256(0x14b36066f483c881dbeda713a2f028597c0d56eea71c2c7e9cf1bf9414bdfafb), uint256(0x0eb92f60db88dc730d09ba6fda75b67f74949a2d2e6608ae07049096a8ea2544));
        vk.gammaABC[21] = Pairing.G1Point(uint256(0x13b87f5d4ba3dc6bfdc6302a45e1bcdf346cf35f3da558599340dde8475ebd0e), uint256(0x25ddd67ec22e8c3e39620e6183e10f4dd0edc82ea9d5e4edaaa69563a218502e));
        vk.gammaABC[22] = Pairing.G1Point(uint256(0x2de068cc97683feb01e58ee6a275ec497b8efcf10d65ce00ff2a5ec5ac787dbc), uint256(0x26c30e60ce44a65b47516cbd55a0ff371d082dab975349651e4d5a19e11b0eb9));
        vk.gammaABC[23] = Pairing.G1Point(uint256(0x2d6489b214b222f6850a51fe99f2a696137f42058e96181abea6818854cd5696), uint256(0x09f46e9580ea5fc0bd0f91e5f0d2dd458ad41f00ee674ac1afbb69a5c3a1687f));
        vk.gammaABC[24] = Pairing.G1Point(uint256(0x1f0401e7a9a5ee66d7fc91704e7340feff821f1447ca23ddafe37a4ddebcfe0e), uint256(0x27a91b7ad861c81f66dfe429ae7990509ada9306a72c162374ae3939c2c50c50));
        vk.gammaABC[25] = Pairing.G1Point(uint256(0x23710013291407d74e605a055991eb7dd89b5ad6d2d082af4b380adff8cb2f14), uint256(0x246a22270dec5a08f2fe8faab620a30531a937945a31db2fbaa779e837ace7bf));
        vk.gammaABC[26] = Pairing.G1Point(uint256(0x0f4ebd9dc72153536a758712f12ebf22beb8492c1faaf4e12392a2c034640413), uint256(0x2db2e7c864e72c62448b0dc2bdfa33bb72d142d5a8b5bbc598c046ea8e545040));
        vk.gammaABC[27] = Pairing.G1Point(uint256(0x092939ceeb4eb975a5357273380882b6c03b90e9398cc513293df1a8f9343618), uint256(0x0cf8868558a1374681a5a7d9df52d76639cde91c5bd6fe0b2315ef6872acb869));
        vk.gammaABC[28] = Pairing.G1Point(uint256(0x030ecce965171577cc0a7bb12afc82cb076fd52ade90c52467f6a316171a202b), uint256(0x14eed2bd3544b7070a7d07bbb738e8f0eba91d49e35250d3bc0f7919b674d913));
        vk.gammaABC[29] = Pairing.G1Point(uint256(0x0572f9b1ffbc6457deba5b28de4a8c5d7948e6c63772a228acc8a8b02118bba6), uint256(0x1f9af15f9cfbb59bc58c6cbc07c1d474762b16da6b8e5f8691381c4ec8054d78));
        vk.gammaABC[30] = Pairing.G1Point(uint256(0x0d589f7c7d462d7c3bbd1e3d181ac7d679775d0f3fa32a1678974c52285a10cf), uint256(0x07b922483132a30a23ee44487548e435a7430cd1d914d8081e31e3c09fbcce7d));
        vk.gammaABC[31] = Pairing.G1Point(uint256(0x149733f67e0015face256038ca413312e935058a58df76f3322ec86deeeb72a6), uint256(0x1315049e2ea5aca332a0732a6d266684ef0be56d8f99e1d19df6304b30d5807f));
        vk.gammaABC[32] = Pairing.G1Point(uint256(0x15656ab57cf9a3d643167ea04be04490b2d919cfe3262d8af2680c0d4d298113), uint256(0x051df09eadf5ca8640438188f7bd46060666ac5aa03ef695e141b0daf59c8634));
        vk.gammaABC[33] = Pairing.G1Point(uint256(0x0e3134702881c7ea5833eec3d7566a45ea76e642518c8b7e92f0f6e57ae230b8), uint256(0x1ee9af5b46436f523e5c19e6553e9c3902c9c4c45b16f49da2db00a2e05e8ecb));
        vk.gammaABC[34] = Pairing.G1Point(uint256(0x02f5ca7619984e01efd2a65d146c1af0f89560c811f20fadc852a733c83e252e), uint256(0x264d5c8b4671c182d418f40f5656dc34eefc5850ef6e5d38a9dc5947b6393a40));
        vk.gammaABC[35] = Pairing.G1Point(uint256(0x247033627bb8ea9c49a6c2c7967ad180737716a30e0ea40a8a088291bbf4a3c0), uint256(0x080f6d61d730402c4b7f073b3d1c37358509886991d25c06641ff7746d85da19));
        vk.gammaABC[36] = Pairing.G1Point(uint256(0x0319b962ee4c0b4fc04db903e6e9ba1ce46d5956b13b5c967c845dc15baf09db), uint256(0x240bc05301df72fa98c9770eb8483cab8235cdaceb0247f7808b5d3cba20431a));
        vk.gammaABC[37] = Pairing.G1Point(uint256(0x00c437ed3c67020617db90e67591c10551fd4c9ce882b4e3560624f9f422502e), uint256(0x2add2a89afd29bccb56a48d3cbbb2360d2f430f6c43182f0c4973fef5a5878b1));
        vk.gammaABC[38] = Pairing.G1Point(uint256(0x16ae663ec71e80e4d106dcb5fc24f079e356f49aa9e1346b8f67a9d39bd5fbca), uint256(0x0dafd12060a1cc82fd84dfe617bfb89d17b42851812129c90b1e6aa2adbae787));
    }
    function verify(uint[] memory input, Proof memory proof) internal returns (uint) {
        VerifyingKey memory vk = verifyingKey();
        require(input.length + 1 == vk.gammaABC.length);
        // Compute the linear combination vk_x
        Pairing.G1Point memory vk_x = Pairing.G1Point(0, 0);
        for (uint i = 0; i < input.length; i++)
            vk_x = Pairing.addition(vk_x, Pairing.scalar_mul(vk.gammaABC[i + 1], input[i]));
        vk_x = Pairing.addition(vk_x, vk.gammaABC[0]);
        if(!Pairing.pairingProd4(
             proof.A, proof.B,
             Pairing.negate(vk_x), vk.gamma,
             Pairing.negate(proof.C), vk.delta,
             Pairing.negate(vk.a), vk.b)) return 1;
        return 0;
    }
    event Verified(string s);
    function verifyTx(
            uint[2] memory a,
            uint[2][2] memory b,
            uint[2] memory c,
            uint[38] memory input
        ) public returns (bool r) {
        Proof memory proof;
        proof.A = Pairing.G1Point(a[0], a[1]);
        proof.B = Pairing.G2Point([b[0][0], b[0][1]], [b[1][0], b[1][1]]);
        proof.C = Pairing.G1Point(c[0], c[1]);
        uint[] memory inputValues = new uint[](input.length);
        for(uint i = 0; i < input.length; i++){
            inputValues[i] = input[i];
        }
        if (verify(inputValues, proof) == 0) {
            emit Verified("Transaction successfully verified.");
            return true;
        } else {
            return false;
        }
    }
}
