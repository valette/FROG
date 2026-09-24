# FROG Codebase Analysis: Bugs, Vulnerabilities & Enhancements

This document summarizes the comprehensive audit of the **FROG** (Fast Registration Of image Groups) codebase, identifying critical bugs, algorithmic issues, command-line parsing vulnerabilities, workflow/scripting flaws, build system gaps, and recommended enhancements, along with resolution status and commit tracking.

---

## Audit Status Dashboard

| ID | Category | Summary | Status | Commit |
|:---|:---|:---|:---|:---|
| **BUG-01** | Memory & Pointer | Pass-by-value prevents points coordinate transformation in `match` | **FIXED** | `074bcef` |
| **BUG-02** | Memory & Pointer | Use-after-free and memory leak in `readFrogJSON` | **FIXED** | `25b31a5` |
| **BUG-03** | Memory & Pointer | `AverageVolumes` accumulates on uninitialized memory for standard deviation | **FIXED** | `663109c` |
| **BUG-04** | Numerical | Floating-point cancellation produces `NaN` in standard deviation | **FIXED** | `663109c` |
| **BUG-05** | Memory & Pointer | Uninitialized pointer dereference on unrecognized output extension | **FIXED** | `80cb6f3` |
| **BUG-06** | Portability | Non-standard Variable Length Arrays (VLAs) in C++ | **FIXED** | `f4e5dea` |
| **BUG-07** | Bounds / Buffer | Out-of-bounds memory writes in `saveErrorMaps` | **FIXED** | `a69a9e7` |
| **BUG-08** | Python Scope | Missing module scope `NameError` in `tools/transformIO.py` | **FIXED** | `39706e6` |
| **BUG-09** | Python Syntax | Python 2 print syntax error in `tools/checkerBoard.py` | **FIXED** | `155ac70` |
| **BUG-10** | Resource Leak | Pervasive VTK and raw pointer leaks across tools | **FIXED** | `e785c87` |
| **ALG-01** | Algorithmic | Inverted fixed image indexing in `setupLinearTransforms` | **FIXED** | `bbb9ba0` |
| **ALG-02** | Numerical | Unprotected denominator in EM distribution estimation | **FIXED** | `5e852c0` |
| **ALG-03** | Robustness | Zero samples causes crash / undefined behavior in `Stats` | **FIXED** | `e332b4f` |
| **ALG-04** | Statistical | Severe sampling bias in `Stats::addSample` | **FIXED** | `ddfc269` |
| **ALG-05** | Algorithmic | Infinite loop risk and coordinate mismatch in `RANSACBatch` | **FIXED** | `b1bf662` |
| **CLI-01** | CLI Safety | Universal out-of-bounds `argv` access across all CLI tools | **FIXED** | `e972c68` |
| **CLI-02** | CLI Bug | `PointsTransform` ignores `-o outputFileName` | **FIXED** | `f1eb113` |
| **CLI-03** | Portability | Missing `<iomanip>` header in `CheckDiffeomorphism` | **FIXED** | `a545f7a` |
| **SCR-01** | Scripting | `chdir(outputDirectory)` invalidates all relative input paths | **FIXED** | `920d38d` |
| **SCR-02** | Scripting | Shared global filename `RAS.nii.gz` overwrites inputs in `FROG.py` | **FIXED** | `eef6637` |
| **SCR-03** | Scripting | Dead, broken copy-paste code in `tools/register2volumes.py` | **FIXED** | `7247305` |
| **SCR-04** | Scripting | Duplicate command line flags in `tools/register2volumes.py` | **FIXED** | `7247305` |
| **SCR-05** | Scripting | Outdated `.tfm` format in `checkDiffeomorphism.sh` | **FIXED** | `a74c658` |
| **SCR-06** | Security | Unquoted command execution and shell injection vulnerability | **FIXED** | `d5570cf` |
| **BLD-01** | Build | Overly restrictive `cmake_minimum_required(VERSION 3.30)` | **FIXED** | `9966589` |
| **BLD-02** | Build | Reference to non-existent `ConvertFrogTransform.cxx` | **RECONCILED** | Audit resolution |
| **BLD-03** | Build / ODR | Non-`inline` functions in `tools/transformIO.h` | **FIXED** | `8af18ca` |
| **BLD-04** | Build | Hardcoded `-mavx` and deprecated `stdc++fs` | **FIXED** | `a52e22d` |
| **BLD-05** | CI / Testing | CI configures testing but defines no tests | **FIXED** | `1b183fe` |
| **ENH-04** | Testing | Add unit test suite (Stats, transformIO) and CI ctest step | **FIXED** | `1b183fe` |
| **ENH-05** | Compliance | Align `package.json` license with `LICENCE.txt` (CeCILL-B) | **FIXED** | Current commit |

---

## 1. Critical Bugs & Memory Errors

### [BUG-01] Pass-by-value prevents points coordinate transformation in `match`
- **Location**: `match/match.cpp`
- **Status**: **FIXED** (`074bcef`)
- **Problem**:
  ```cpp
  for ( auto point : *points ) {
      if ( transform ) {
          transform->TransformPoint( point.coordinates, point.transformedCoordinates );
      } else {
          for ( int i = 0; i < 3; i++ )
              point.transformedCoordinates[ i ] = point.coordinates[ i ];
      }
  }
  ```
  `for ( auto point : *points )` copied `Point` by value. Modifying `point.transformedCoordinates` only mutated the local temporary.
- **Fix**: Use reference: `for ( auto &point : *points )`.

---

### [BUG-02] Use-after-free and memory leak in `readFrogJSON`
- **Location**: `tools/transformIO.h`
- **Status**: **FIXED** (`25b31a5`)
- **Problem**: Accessing `reader->GetQFormMatrix()` after `reader->Delete()`.
- **Fix**: Reordered access and adopted `vtkSmartPointer<vtkNIFTIImageReader>`.

---

### [BUG-03] `AverageVolumes` accumulates on uninitialized memory for standard deviation
- **Location**: `tools/AverageVolumes.cxx`
- **Status**: **FIXED** (`663109c`)
- **Problem**: `ptrOut2` was initialized to `average->GetScalarPointer()` twice, leaving `stdev` memory uninitialized before accumulation.
- **Fix**: Initialized `ptrOut2` to `stdev->GetScalarPointer()`.

---

### [BUG-04] Floating-point cancellation produces `NaN` in standard deviation
- **Location**: `tools/AverageVolumes.cxx`
- **Status**: **FIXED** (`663109c`)
- **Problem**: In constant regions where $E[X^2] \approx (E[X])^2$, rounding errors could make difference slightly negative, resulting in `sqrt(negative) == NaN`.
- **Fix**: Clamped variance using `std::max(0.0f, diff)` prior to `sqrt`.

---

### [BUG-05] Uninitialized pointer dereference on unrecognized output extension
- **Location**: `tools/VolumeTransform.cxx`
- **Status**: **FIXED** (`80cb6f3`)
- **Problem**: If an unsupported extension was passed, `writer` remained an uninitialized pointer and was dereferenced at `writer->SetFileName()`.
- **Fix**: Initialized `writer = nullptr` and added validation exiting cleanly with an error message.

---

### [BUG-06] Non-standard Variable Length Arrays (VLAs) in C++
- **Location**: `registration/imageGroup.cxx`
- **Status**: **FIXED** (`f4e5dea`)
- **Problem**: C99 VLAs (`float *coeffs[nImages];`, `float anchors[...][3];`, `char name[nameLength + 1];`) break MSVC and risk stack overflows.
- **Fix**: Replaced with `std::vector<float*>`, `std::vector<std::array<float, 3>>`, and `std::string`.

---

### [BUG-07] Out-of-bounds memory writes in `saveErrorMaps`
- **Location**: `registration/imageGroup.cxx`
- **Status**: **FIXED** (`a69a9e7`)
- **Problem**: Point coordinates shifted outside the volume bounds yielded out-of-range voxel indices, corrupting heap memory.
- **Fix**: Added boundary checks verifying `intCoord[k] >= 0 && intCoord[k] < dims[k]` for all dimensions before indexing.

---

### [BUG-08] Missing module scope `NameError` in `tools/transformIO.py`
- **Location**: `tools/transformIO.py`
- **Status**: **FIXED** (`39706e6`)
- **Problem**: Unqualified `vtkMatrix4x4()` when module was imported as `import vtk`.
- **Fix**: Changed to `vtk.vtkMatrix4x4()`.

---

### [BUG-09] Python 2 print syntax error in `tools/checkerBoard.py`
- **Location**: `tools/checkerBoard.py`
- **Status**: **FIXED** (`155ac70`)
- **Problem**: Missing parentheses `print pattern` caused Python 3 `SyntaxError`.
- **Fix**: Converted to `print(pattern)`.

---

### [BUG-10] Pervasive VTK and raw pointer leaks across tools
- **Location**: `tools/MeshTransform.cxx`, `tools/PointsTransform.cxx`, `tools/CheckDiffeomorphism.cxx`, `tools/DummyVolumeGenerator.cxx`, `tools/VolumeTransform.cxx`, `match/match.cpp`
- **Status**: **FIXED** (`e785c87`)
- **Problem**: VTK objects allocated with `::New()` and raw heap buffers were never released.
- **Fix**: Added proper `Delete()` calls and converted objects to RAII / `vtkSmartPointer`. Closed opened `FILE*` handles.

---

## 2. Algorithmic & Statistical Defects

### [ALG-01] Inverted fixed image indexing in `setupLinearTransforms`
- **Location**: `registration/imageGroup.cxx`
- **Status**: **FIXED** (`bbb9ba0`)
- **Problem**: Fixed images are stored at indices `0 .. numberOfFixedImages - 1`, but anchor averaging checked `i < images.size() - numberOfFixedImages`.
- **Fix**: Checked `i < this->numberOfFixedImages` correctly and safeguarded against zero division.

---

### [ALG-02] Unprotected denominator in EM distribution estimation
- **Location**: `registration/stats.cxx`
- **Status**: **FIXED** (`5e852c0`)
- **Problem**: Denominator `sum4` was not protected with `epsilon`, producing `NaN` or `Inf` if all points were inliers.
- **Fix**: Protected denominator with `sum4 = std::max(sum4, epsilon)`.

---

### [ALG-03] Zero samples causes crash / undefined behavior in `Stats`
- **Location**: `registration/stats.cxx`
- **Status**: **FIXED** (`e332b4f`)
- **Problem**: Calling `estimateDistribution()`, `getHistogram()`, or `displayParameters()` on empty stats dereferenced end iterators or divided by zero.
- **Fix**: Added `if (this->size == 0) return;` safeguards across `Stats` methods.

---

### [ALG-04] Severe sampling bias in `Stats::addSample`
- **Location**: `registration/stats.h`
- **Status**: **FIXED** (`ddfc269`)
- **Problem**: Sample collection stopped accepting new samples once `maxSize` was reached, biasing distribution fitting toward the first points encountered.
- **Fix**: Implemented Reservoir Sampling (Algorithm R) with uniform probability across all processed points.

---

### [ALG-05] Potential infinite loop and coordinate space mismatch in `RANSACBatch`
- **Location**: `registration/imageGroup.cxx`
- **Status**: **FIXED** (`b1bf662`)
- **Problem**: Loop `while (true)` spun indefinitely if points had no links; modulo zero if `nPoints == 0`; coordinate space discrepancy between target landmarks (`xyz`) and consensus evaluation (`xyz2`).
- **Fix**: Validated candidate pool, guarded iteration count, and unified landmark coordinate spaces.

---

## 3. Command-Line & Argument Parsing Vulnerabilities

### [CLI-01] Universal out-of-bounds `argv` access across all CLI tools
- **Locations**: `registration/frog.cxx`, `match/match.cpp`, `tools/VolumeTransform.cxx`, `tools/MeshTransform.cxx`, `tools/PointsTransform.cxx`
- **Status**: **FIXED** (`e972c68`)
- **Problem**: Unchecked `argv[argumentsIndex + 1]` triggered segfault when options were passed without values.
- **Fix**: Added bounds checking `argumentsIndex + 1 >= argc` with helpful error diagnostics.

---

### [CLI-02] `PointsTransform` ignores `-o outputFileName`
- **Location**: `tools/PointsTransform.cxx`
- **Status**: **FIXED** (`f1eb113`)
- **Problem**: Parsed `-o` flag but never wrote coordinates to the specified output file.
- **Fix**: Implemented file writing for transformed coordinates when `-o` is supplied.

---

### [CLI-03] Missing `<iomanip>` header in `CheckDiffeomorphism`
- **Location**: `tools/CheckDiffeomorphism.cxx`
- **Status**: **FIXED** (`a545f7a`)
- **Problem**: Used `std::setprecision` without direct include.
- **Fix**: Added `#include <iomanip>`.

---

## 4. Scripting & Workflow Defects (Python & Shell)

### [SCR-01] `chdir(outputDirectory)` invalidates all relative input paths
- **Location**: `FROG.py`, `run.sh`, `transform.sh`
- **Status**: **FIXED** (`920d38d`)
- **Problem**: Relative paths in volume lists broke when script changed working directory.
- **Fix**: Converted all paths to absolute paths before calling `os.chdir()` or `cd`.

---

### [SCR-02] Shared global filename `RAS.nii.gz` overwrites inputs in `FROG.py`
- **Location**: `FROG.py`
- **Status**: **FIXED** (`eef6637`)
- **Problem**: Flipped volumes wrote to identical filename `RAS.nii.gz`, repeatedly overwriting each other.
- **Fix**: Generated unique per-image filenames `f"RAS_{i}.nii.gz"`.

---

### [SCR-03] Dead, broken copy-paste code in `tools/register2volumes.py`
- **Location**: `tools/register2volumes.py`
- **Status**: **FIXED** (`7247305`)
- **Problem**: Function `computeAverageImage` referenced nonexistent variables and crashed.
- **Fix**: Removed dead copy-pasted function.

---

### [SCR-04] Duplicate command line flags in `tools/register2volumes.py`
- **Location**: `tools/register2volumes.py`
- **Status**: **FIXED** (`7247305`)
- **Problem**: `-s` and `-p` flags were appended twice to argument list.
- **Fix**: Cleaned argument builder and used `sys.executable`.

---

### [SCR-05] Outdated `.tfm` format in `checkDiffeomorphism.sh`
- **Location**: `checkDiffeomorphism.sh`
- **Status**: **FIXED** (`a74c658`)
- **Problem**: Script pointed to obsolete `tfm/$CUR_IT.tfm` instead of JSON transforms.
- **Fix**: Updated path to `transforms/$CUR_IT.json`.

---

### [SCR-06] Unquoted command execution and shell injection vulnerability
- **Location**: `FROG.py`, `tools/register.py`, `tools/register2volumes.py`, `run.sh`
- **Status**: **FIXED** (`d5570cf`)
- **Problem**: Vulnerable string concatenations executed via `os.system()` and unquoted bash arguments.
- **Fix**: Migrated to `subprocess.run(..., check=True)` with argument lists; quoted bash variables.

---

## 5. Build System, CMake & CI/CD Gaps

### [BLD-01] Overly restrictive `cmake_minimum_required(VERSION 3.30)`
- **Location**: `CMakeLists.txt`
- **Status**: **FIXED** (`9966589`)
- **Fix**: Lowered requirement to `cmake_minimum_required(VERSION 3.20)` for LTS compatibility.

---

### [BLD-02] Audit Reconciliation: `tools/ConvertFrogTransform.cxx`
- **Location**: `tools/CMakeLists.txt`
- **Status**: **RECONCILED** (Audit Finding Resolved)
- **Note**: The audit originally noted `tools/ConvertFrogTransform.cxx` missing from CMake. Verification showed that this file never existed in the repository or git history; it was an audit hallucination. No phantom file or CMake entry was added.

---

### [BLD-03] Non-`inline` functions in `tools/transformIO.h` (ODR / Linker risk)
- **Location**: `tools/transformIO.h`
- **Status**: **FIXED** (`8af18ca`)
- **Fix**: Marked free helper functions with `inline` to prevent ODR duplicate symbol linker errors.

---

### [BLD-04] Hardcoded `-mavx` and deprecated `stdc++fs`
- **Location**: `match/CMakeLists.txt`, `registration/CMakeLists.txt`
- **Status**: **FIXED** (`a52e22d`)
- **Fix**: Made AVX flags portable across MSVC and GCC/Clang (`/arch:AVX` vs `-mavx`). Linked `stdc++fs` conditionally only for GCC < 9.0.

---

### [BLD-05] CI configures testing but defines no tests and runs no tests
- **Location**: `.github/workflows/ci.yml`, `CMakeLists.txt`, `tests/`
- **Status**: **FIXED** (`1b183fe`)
- **Fix**: Added CTest test suite with unit tests in `tests/`, enabled testing in CMake, and updated GitHub Actions CI with `ctest --output-on-failure` and `actions/checkout@v4`.

---

## 6. Enhancements & Modernization Roadmap

- **[ENH-01] Architecture & Robustness**: Safe CLI parsing implemented across tools; VTK memory management cleaned up; boundary checks added.
- **[ENH-02] Performance & Scalability**: Reservoir sampling added to `Stats` (Algorithm R); SIMD flags made portable.
- **[ENH-03] Python Tools Modernization**: Full migration to `subprocess.run`; safe argument passing and path handling implemented.
- **[ENH-04] Testing & Quality Assurance**: Unit tests for `Stats` and `transformIO` added to `tests/` and integrated into CTest and CI.
- **[ENH-05] Documentation & License Consistency**: `package.json` license aligned to `CECILL-B` matching `LICENCE.txt` and `Readme.md`.
