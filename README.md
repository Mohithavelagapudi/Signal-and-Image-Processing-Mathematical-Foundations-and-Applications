# Signal and Image Processing: Mathematical Foundations and Applications

### Overview

Signal and image processing are foundational pillars in modern data science, AI, and engineering. This repository presents a curated set of MATLAB codes and experiments that demonstrate the mathematical underpinnings and practical applications of core techniques such as Fourier analysis, convolution, clustering, and compression. Each script is designed to bridge theory and practice, providing both a learning resource and a toolkit for real-world problem solving in domains like Computer Vision, Audio Analytics, and Pattern Recognition.

---
# 🧠 Mathematical Foundations of Signal & Image Processing

---

### 📐 Mathematical Foundations

| **Technique** | **Mathematical Principle** | **Key Equations / Concepts** |
|---------------|-----------------------------|-------------------------------|
| **Discrete Fourier Transform (DFT)** | Decomposition of signals into frequency components | $$ X_k = \sum_{n=0}^{N-1} x_n \, e^{-j 2 \pi kn / N} $$ |
| **Convolution** | Linear system response, time/frequency domain duality | $$ (f * g)(t) = \int f(\tau) \, g(t - \tau) \, d\tau $$ |
| **Image DFT & Processing** | 2D Fourier Transform for spatial frequency analysis | $$ F(u,v) = \sum_{x} \sum_{y} f(x,y) \, e^{-j 2\pi(ux/M + vy/N)} $$ |
| **K-Means Clustering** | Unsupervised learning — minimization of within-cluster variance | $$ \min \sum_{i=1}^{k} \sum_{x \in C_i} \| x - \mu_i \|^2 $$ |
| **Kronecker & Hadamard Products** | Matrix operations for multi-dimensional data manipulation | $$ \text{Kronecker: } (A \otimes B), \quad \text{Hadamard: } (A \circ B) $$ |
| **Fourier Series** | Periodic signal decomposition into sines and cosines | $$ f(t) = a_0 + \sum_{n=1}^{N} [a_n \cos(n\omega_0 t) + b_n \sin(n\omega_0 t)] $$ |

---

### 🧩 Project Structure & Highlights

#### **1️⃣ Fourier Analysis & Signal Processing**

| **File** | **Description** |
|-----------|----------------|
| `signalanalysis_discrete.m` | Manual DFT computation and spectral visualization. |
| `fouriertransform_signalprocessing.m` | Implements frequency-domain filtering and reconstruction. |
| `derivative_ft.m` | Numerical differentiation using Fourier derivative properties. |
| `convolution.m` | Demonstrates time-domain convolution using Fourier equivalence. |

---

#### **2️⃣ Image Processing & Compression**

| **File** | **Description** |
|-----------|----------------|
| `image_processing.m` | 2D DFT computation, spatial filtering, and visualization of magnitude spectra. |
| `imagecompression.m` | Image compression via thresholding of frequency coefficients (energy compaction). |

---

### ⚙️ Mathematical Insights

#### **Fourier Differentiation Property**

The derivative of a signal in the time domain corresponds to multiplication by \( j\omega \) in the frequency domain:

$$
\frac{d f(t)}{d t} \; \xleftrightarrow{\mathcal{F}} \; j \omega F(\omega)
$$

This property enables **numerical differentiation** without finite difference approximations.

---

#### **Convolution Theorem**

Convolution in the time domain is equivalent to multiplication in the frequency domain:

$$
f(t) * g(t) \; \xleftrightarrow{\mathcal{F}} \; F(\omega) \cdot G(\omega)
$$

This allows **filtering and signal transformation** to be efficiently performed in the frequency domain.

---

#### **Energy Compaction (Image Compression)**

Most of an image’s energy is concentrated in **low-frequency components**.  
By discarding high-frequency coefficients, we achieve **lossy compression** with minimal perceptual degradation:

$$
E_{\text{low}} = \sum_{(u,v) \in \text{low-freq}} |F(u,v)|^2
$$

Thus, retaining a small subset of low-frequency terms preserves most visual information.

---
### 🔑 Getting Started

- Requirements: MATLAB (or GNU Octave for most scripts), sample images/audio files as referenced in the code.
- Usage: Run individual .m files to explore each concept. Modify parameters (e.g., compression ratio, number of clusters) to observe effects.
- Visualization: Most scripts include plotting for intuitive understanding.

---
### 🧠 Acknowledgements

This repository is the result of hands-on experimentation and learning in the field of signal and image processing, with inspiration from academic coursework and open-source contributions.

📚 *Developed as part of the Signal & Image Processing Lab Exercises (MATLAB).*
