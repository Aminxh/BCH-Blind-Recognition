# BCH-Blind-Recognition

A **C++ implementation of blind recognition for BCH codes**. Automatically estimates BCH code parameters such as codeword length, number of information bits, error-correcting capability, and primitive polynomial from encoded data without prior knowledge of the code structure.

---

## Features

- Blind recognition of BCH codes from input binary data.
- Supports both **extended** and **shortened** BCH codes.
- Automatically identifies:
  - Codeword length `n`
  - Number of information bits `k`
  - Error-correcting capability `t`
  - Primitive polynomial used
- Threshold-based detection to improve reliability.
- Uses **Armadillo** for efficient matrix operations.

---

## Requirements

- **C++17** or higher  
- [Armadillo](http://arma.sourceforge.net/) linear algebra library  
- Standard C++ library: `<iostream>`  

Optional: Precomputed Galois Field tables included (`GF_TABLE1.h`, `GF_TABLE2.h`, `GF_TABLE_PRIM_POLY.h`).

---

## File Structure

```bash
├── include/                   # Header files
│   ├── BCH_Blind.h            # Header file with class and structure definitions
│   ├── GF_TABLE1.h            # Galois field precomputed table 1
│   ├── GF_TABLE2.h            # Galois field precomputed table 2
│   └── GF_TABLE_PRIM_POLY.h   # Galois field primitive polynomials
│
├── src/                       # Source code files
│   ├── BCH_Blind.cpp          # Implementation of BCH blind recognition
│   └── main.cpp               # example usage

```


## Algorithm Overview

1. Input data is reshaped according to candidate codeword lengths.

2. The algorithm iterates over all possible primitive polynomials for the given codeword length.

3. Uses Galois Field arithmetic to compute syndrome-like patterns.

4. Detects repeated error patterns to estimate the number of correctable errors (t).

5. Returns estimated BCH parameters if recognition is successful.

## Notes

Only supports specific BCH code lengths (n ≤ 255) as per myBCHNumerr() table.

Input data must be in uchar_rowvec format (0s and 1s).

Threshold controls the strictness of recognition; typical values: 0.6–0.8.

Only the first 1000 messages are used for recognition for efficiency.

## References

Lin, Shu, and Daniel J. Costello. Error Control Coding: Fundamentals and Applications. 2nd Edition.

[BCH code – Wikipedia](https://en.wikipedia.org/wiki/BCH_code)
