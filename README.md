# Optics.jl
Analyze ray traces in a simple syntax with [Julia Lang](https://julialang.org)

⚠️**Under development.** Almost functionality to be included in this package still not implemented.

<img src="https://github.com/user-attachments/assets/bb928571-2196-4ea8-9f22-832648d90d93">

## Features
### Rays
* [x] **Single Ray**
* [x] **Point Source**

### Spherical Lenses
* [x] **Plano-Convex Lens**
<img src="https://github.com/user-attachments/assets/f234fffc-4a38-4f11-b046-bec8ebc3ee02" width="50%">

* [x] **Plano-Convex Lens**
<img src="https://github.com/user-attachments/assets/2189526a-42d9-4144-9ba8-29975c0d9425" width="50%">

* [x] **Bi-Convex Lens**
<img src="https://github.com/user-attachments/assets/288f930b-7697-48a0-8f4b-52f1515bb724" width="50%">

* [ ] **Bi-Concave Lens**
* [ ] **Positive Meniscus Lens**
* [ ] **Negative Meniscus Lens**

### Other Lenses
* [x] **Aspheric Convex Lenss**
<img src="https://github.com/user-attachments/assets/38919070-402c-49fb-ab62-632f295c0341" width="50%">

  ```julia
  # Argument orderes =====================================================================================
  # diameter, thickness, left aspheric curve, right aspheric curve, refractive index, a flag for mirrored
  # ======================================================================================================
  Optics.AsphericConvexLens(25.4, 14.0,
    # diameter, coefficients, conic constant, radius of curvature
    Optics.AsphericCurve(25.4, [8.6821674e-05, 6.3760123e-08, 2.4073084e-09, -1.7189021e-011], -0.9991715, 8.818197),
    Optics.AsphericCurve(25.4, [0.0], 0.0, -69.99948),
    1.52,
    is_mirrored=true
  )
  ```

* [ ] **Fresnel Lens**


### Mirrors
* [ ] **Concave Mirror**

### Miscellaneous Elements
* [x] **Wall**
