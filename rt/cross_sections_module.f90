! cross_sections_module.f90
module cross_sections_module
  implicit none

  private
  public::initialize_cross_sections, getCrosssection_rtz

  type cross_sections
      real(KIND=8)   :: E_th(27,27)
      real(KIND=8)   :: E_max(27,27)
      real(KIND=8)   :: E_0(27,27)
      real(KIND=8)   :: sig_0(27,27)
      real(KIND=8)   :: y_a(27,27)
      real(KIND=8)   :: P(27,27)
      real(KIND=8)   :: y_w(27,27)
      real(KIND=8)   :: y_0(27,27)
      real(KIND=8)   :: y_1(27,27)
  end type cross_sections

  type(cross_sections) :: verner_cross_sections

CONTAINS

SUBROUTINE initialize_cross_sections()
  ! Photoionization cross sections from
  ! https://articles.adsabs.harvard.edu/pdf/1996ApJ...465..487V
  implicit none

  verner_cross_sections%E_th = 0.d0
  verner_cross_sections%E_max = 0.d0
  verner_cross_sections%E_0 = 0.d0
  verner_cross_sections%sig_0 = 0.d0
  verner_cross_sections%y_a = 0.d0
  verner_cross_sections%P = 0.d0
  verner_cross_sections%y_w = 0.d0
  verner_cross_sections%y_0 = 0.d0
  verner_cross_sections%y_1 = 0.d0
  
  ! Hydrogen
  verner_cross_sections%E_th(1,1:1) = (/ 1.360d1 /)
  verner_cross_sections%E_max(1,1:1) = (/ 5.000d4 /)
  verner_cross_sections%E_0(1,1:1) = (/ 4.298d-1 /)
  verner_cross_sections%sig_0(1,1:1) = (/ 5.475d4 /)
  verner_cross_sections%y_a(1,1:1) = (/ 3.288d1 /)
  verner_cross_sections%P(1,1:1) = (/ 2.963d0 /)
  verner_cross_sections%y_w(1,1:1) = (/ 0.d0 /)
  verner_cross_sections%y_0(1,1:1) = (/ 0.d0 /)
  verner_cross_sections%y_1(1,1:1) = (/ 0.d0 /)

  ! Helium
  verner_cross_sections%E_th(2,1:2) = (/ 2.459d1, 5.442d1 /)
  verner_cross_sections%E_max(2,1:2) = (/ 5.000d4, 5.000d4 /)
  verner_cross_sections%E_0(2,1:2) = (/ 1.361d1, 1.720d0 /)
  verner_cross_sections%sig_0(2,1:2) = (/ 0.492d2, 1.369d4 /)
  verner_cross_sections%y_a(2,1:2) = (/ 1.469d0, 3.288d1 /)
  verner_cross_sections%P(2,1:2) = (/ 3.188d0, 2.963d0 /)
  verner_cross_sections%y_w(2,1:2) = (/ 2.039d0, 0.d0 /)
  verner_cross_sections%y_0(2,1:2) = (/ 4.434d-1, 0.d0 /)
  verner_cross_sections%y_1(2,1:2) = (/ 2.136d0, 0.d0 /)

  ! Carbon
  verner_cross_sections%E_th(6,1:6) = (/ 1.126d1, 2.438d1, 4.789d1, 6.449d1, 3.921d2, 4.900d2 /)
  verner_cross_sections%E_max(6,1:6) = (/ 2.910d2, 3.076d2, 3.289d2, 3.522d2, 5.000d4, 5.000d4 /)
  verner_cross_sections%E_0(6,1:6) = (/ 2.144d0, 4.058d-1, 4.614d0, 3.506d0, 4.624d1, 1.548d1 /)
  verner_cross_sections%sig_0(6,1:6) = (/ 5.027d2, 8.709d0, 1.539d4, 1.068d2, 2.344d2, 1.521d3 /)
  verner_cross_sections%y_a(6,1:6) = (/ 6.216d1, 1.261d2, 1.737d0, 1.436d1, 2.183d1, 3.288d1 /)
  verner_cross_sections%P(6,1:6) = (/ 5.101d0, 8.578d0, 1.593d1, 7.457d0, 2.581d0, 2.963d0 /)
  verner_cross_sections%y_w(6,1:6) = (/ 9.157d-2, 2.093d0, 5.922d0, 0.d0, 0.d0, 0.d0 /)
  verner_cross_sections%y_0(6,1:6) = (/ 1.133d0, 4.929d1, 4.378d-3, 0.d0, 0.d0, 0.d0 /)
  verner_cross_sections%y_1(6,1:6) = (/ 1.607d0, 3.234d0, 2.528d-2, 0.d0, 0.d0, 0.d0 /)

  ! Nitrogen
  verner_cross_sections%E_th(7,1:7) = (/ 1.453d1, 2.960d1, 4.745d1, 7.747d1, 9.789d1, 5.521d2, 6.671d2 /)
  verner_cross_sections%E_max(7,1:7) = (/ 4.048d2, 4.236d2, 4.473d2, 5.753d2, 5.043d2, 5.0d4, 5.0d4 /)
  verner_cross_sections%E_0(7,1:7) = (/ 4.034d0, 6.128d-2, 2.420d-1, 5.494d0, 4.471d0, 6.943d1, 2.108d1 /)
  verner_cross_sections%sig_0(7,1:7) = (/ 8.235d2, 1.944d0, 9.375d-1, 1.690d4, 8.376d1, 1.519d2, 1.117d3 /)
  verner_cross_sections%y_a(7,1:7) = (/ 8.033d1, 8.163d2, 2.788d2, 1.714d0, 3.297d1, 2.627d1, 3.288d1 /)
  verner_cross_sections%P(7,1:7) = (/ 3.928d0, 8.773d0, 9.156d0, 1.706d1, 6.003d0, 2.315d0, 2.963d0 /)
  verner_cross_sections%y_w(7,1:7) = (/ 9.097d-2, 1.043d1, 1.850d0, 7.904d0, 0.d0, 0.d0, 0.d0 /)
  verner_cross_sections%y_0(7,1:7) = (/ 8.598d-1, 4.280d2, 1.877d2, 6.415d-3, 0.d0, 0.d0, 0.d0 /)
  verner_cross_sections%y_1(7,1:7) = (/ 2.325d0, 2.030d1, 3.999d0, 1.937d-2, 0.d0, 0.d0, 0.d0 /)

  ! Oxygen
  verner_cross_sections%E_th(8,1:8) = (/ 1.362d1, 3.512d1, 5.494d1, 7.741d1, 1.139d2, 1.381d2, 7.393d2, 8.714d2 /)
  verner_cross_sections%E_max(8,1:8) = (/ 5.380d2, 5.581d2, 5.840d2, 6.144d2, 6.491d2, 6.837d2, 5.000d4, 5.000d4 /)
  verner_cross_sections%E_0(8,1:8) = (/ 1.240d0, 1.386d0, 1.723d-1, 2.044d-1, 2.854d0, 7.824d0, 8.709d0, 2.754d1 /)
  verner_cross_sections%sig_0(8,1:8) = (/ 1.745d3, 5.967d1, 6.753d2, 8.659d-1, 1.642d4, 6.864d1, 1.329d2, 8.554d2 /)
  verner_cross_sections%y_a(8,1:8) = (/ 3.784d0, 3.175d1, 3.852d2, 4.931d2, 1.792d0, 3.210d1, 2.535d1, 3.288d1 /)
  verner_cross_sections%P(8,1:8) = (/ 1.764d1, 8.943d0, 6.822d0, 8.785d0, 2.647d1, 5.495d0, 2.336d0, 2.963d0 /)
  verner_cross_sections%y_w(8,1:8) = (/ 7.589d-2, 1.934d-2, 1.191d-1, 3.143d0, 2.836d1, 0.d0, 0.d0, 0.d0 /)
  verner_cross_sections%y_0(8,1:8) = (/ 8.698d0, 2.131d1, 3.839d-3, 3.328d2, 3.036d-2, 0d0, 0d0, 0d0 /)
  verner_cross_sections%y_1(8,1:8) = (/ 1.271d-1, 1.503d-2, 4.569d-1, 4.285d1, 5.554d-2, 0d0, 0d0, 0d0 /)

  ! Neon
  verner_cross_sections%E_th(10,1:10) = (/ 2.156d+01, 4.096d+01, 6.346d+01, 9.712d+01, 1.262d+02, 1.579d+02, 2.073d+02, 2.391d+02, 1.196d+03, 1.362d+03 /)
  verner_cross_sections%E_max(10,1:10) = (/ 8.701d+02, 8.831d+02, 9.131d+02, 9.480d+02, 9.873d+02, 1.031d+03, 1.078d+03, 1.125d+03, 5.000d+04, 5.000d+04 /)
  verner_cross_sections%E_0(10,1:10) = (/ 4.870d+00, 1.247d+01, 7.753d-01, 5.566d+00, 1.248d+00, 1.499d+00, 4.888d+00, 1.003d+01, 1.586d+02, 4.304d+01 /)
  verner_cross_sections%sig_0(10,1:10) = (/ 4.287d+03, 1.583d+03, 5.708d+00, 1.685d+03, 2.430d+00, 9.854d-01, 1.198d+04, 5.631d+01, 6.695d+01, 5.475d+02 /)
  verner_cross_sections%y_a(10,1:10) = (/ 5.798d+00, 3.935d+00, 6.725d+01, 6.409d+02, 1.066d+02, 1.350d+02, 1.788d+00, 3.628d+01, 3.352d+01, 3.288d+01 /)
  verner_cross_sections%P(10,1:10) = (/ 8.355d+00, 7.810d+00, 1.005d+01, 3.056d+00, 8.999d+00, 8.836d+00, 2.550d+01, 5.585d+00, 2.002d+00, 2.963d+00 /)
  verner_cross_sections%y_w(10,1:10) = (/ 2.434d-01, 6.558d-02, 4.633d-01, 8.290d-03, 6.855d-01, 1.656d+00, 2.811d+01, 0.000d+00, 0.000d+00, 0.000d+00 /)
  verner_cross_sections%y_0(10,1:10) = (/ 4.236d-02, 1.520d+00, 7.654d+01, 5.149d+00, 9.169d+01, 1.042d+02, 2.536d-02, 0.000d+00, 0.000d+00, 0.000d+00 /)
  verner_cross_sections%y_1(10,1:10) = (/ 5.873d+00, 1.084d-01, 2.023d+00, 6.687d+00, 3.702d-01, 1.435d+00, 4.417d-02, 0.000d+00, 0.000d+00, 0.000d+00 /)

  ! Magnesium
  verner_cross_sections%E_th(12,1:12) = (/ 7.646d0, 1.504d1, 8.014d1, 1.093d2, 1.413d2, 1.865d2, 2.249d2, 2.660d2, 3.282d2, 3.675d2, 1.762d3, 1.963d3 /)
  verner_cross_sections%E_max(12,1:12) = (/ 5.490d1, 6.569d1, 1.317d3, 1.356d3, 1.400d3, 1.449d3, 1.503d3, 1.558d3, 1.618d3, 1.675d3, 5.0d4, 5.0d4 /)
  verner_cross_sections%E_0(12,1:12) = (/ 1.197d1, 8.139d0, 1.086d1, 2.912d1, 9.762d-1, 1.711d0, 3.570d0, 4.884d-1, 3.482d1, 1.452d1, 2.042d2, 6.203d1 /)
  verner_cross_sections%sig_0(12,1:12) = (/ 1.372d8, 3.278d0, 5.377d2, 1.394d3, 1.728d0, 2.185d0, 3.104d0, 6.344d-2, 9.008d2, 4.427d1, 6.140d1, 3.802d2 /)
  verner_cross_sections%y_a(12,1:12) = (/ 2.228d-1, 4.241d7, 9.779d0, 2.895d0, 9.184d1, 9.350d1, 6.060d1, 5.085d2, 1.823d0, 3.826d1, 2.778d1, 3.288d1 /)
  verner_cross_sections%P(12,1:12) = (/ 1.574d1, 3.610d0, 7.117d0, 6.487d0, 1.006d1, 9.202d0, 8.857d0, 9.385d0, 1.444d1, 5.460d0, 2.161d0, 2.963d0 /)
  verner_cross_sections%y_w(12,1:12) = (/ 2.805d-1, 0.d0, 2.604d0, 4.326d-2, 8.090d-1, 6.325d-1, 1.422d0, 6.666d-1, 2.751d0, 0.d0, 0.d0, 0.d0 /)
  verner_cross_sections%y_0(12,1:12) = (/ 0.d0, 0.d0, 4.860d0, 9.402d-1, 1.276d2, 1.007d2, 5.452d1, 5.348d2, 5.444d0, 0.d0, 0.d0, 0.d0 /)
  verner_cross_sections%y_1(12,1:12) = (/ 0.d0, 0.d0, 3.722d0, 1.135d-1, 3.979d0, 1.729d0, 2.078d0, 3.997d-3, 7.918d-2, 0.d0, 0.d0, 0.d0 /)

  ! Silicon
  verner_cross_sections%E_th(14,1:14) = (/ 8.152d0, 1.635d1, 3.349d1, 4.514d1, 1.668d2, 2.051d2, 2.465d2, 3.032d2, 3.511d2, 4.014d2, 4.761d2, 5.235d2, 2.438d3, 2.673d3 /)
  verner_cross_sections%E_max(14,1:14) = (/ 1.060d2, 1.186d2, 1.311d2, 1.466d2, 1.887d3, 1.946d3, 2.001d3, 2.058d3, 2.125d3, 2.194d3, 2.268d3, 2.336d3, 5.0d4, 5.0d4 /)
  verner_cross_sections%E_0(14,1:14) = (/ 2.317d1, 2.556d0, 1.659d-1, 1.288d1, 7.761d-1, 6.305d1, 3.277d-1, 7.655d-1, 3.343d-1, 8.787d-1, 1.205d1, 3.560d1, 2.752d2, 8.447d1 /)
  verner_cross_sections%sig_0(14,1:14) = (/ 2.506d1, 4.140d0, 5.790d-4, 6.083d0, 8.863d-1, 7.293d1, 6.680d-2, 3.477d-1, 1.465d-1, 1.950d-1, 1.992d4, 2.539d1, 4.754d1, 2.793d2 /)
  verner_cross_sections%y_a(14,1:14) = (/ 2.057d1, 1.337d1, 1.474d2, 1.356d6, 1.541d2, 1.558d2, 4.132d1, 3.733d2, 1.404d3, 7.461d2, 1.582d0, 3.307d1, 2.848d1, 3.288d1 /)
  verner_cross_sections%P(14,1:14) = (/ 3.546d0, 1.191d1, 1.336d1, 3.353d0, 9.980d0, 2.400d0, 1.606d1, 8.986d0, 8.503d0, 8.302d0, 2.425d1, 4.728d0, 2.135d0, 2.963d0 /)
  verner_cross_sections%y_w(14,1:14) = (/ 2.837d-1, 1.570d0, 8.626d-1, 0.d0, 1.303d0, 2.989d-3, 3.280d0, 1.476d-3, 1.646d0, 4.489d-1, 2.392d1, 0.d0, 0.d0, 0.d0 /)
  verner_cross_sections%y_0(14,1:14) = (/ 1.672d-5, 6.634d0, 9.613d1, 0.d0, 2.009d2, 1.115d0, 1.149d-2, 3.850d2, 1.036d3, 4.528d2, 1.990d-2, 0.d0, 0.d0, 0.d0 /)
  verner_cross_sections%y_1(14,1:14) = (/ 4.207d-1, 1.272d-1, 6.442d-1, 0.d0, 4.537d0, 8.051d-2, 6.396d-1, 8.999d-2, 2.936d-1, 1.015d0, 1.007d-2, 0.d0, 0.d0, 0.d0 /)

  ! Sulfur
  verner_cross_sections%E_th(16,1:16) = (/ 10.36d0, 23.33d0, 34.83d0, 47.31d0, 72.68d0, &
                                            88.05d0, 280.9d0, 328.2d0, 379.1d0, 447.1d0, &
                                            504.8d0, 564.7d0, 651.7d0, 707.2d0, 3224.d0, &
                                            3494.d0 /)
  verner_cross_sections%E_max(16,1:16) = (/ 170.d0, 184.6d0, 199.5d0, 216.4d0, 235.d0, &
                                            255.7d0, 2569.d0, 2641.d0, 2705.d0, 2782.d0, &
                                            2859.d0, 2941.d0, 3029.d0, 3107.d0, 50000.d0, &
                                            50000.d0 /)
  verner_cross_sections%E_0(16,1:16) = (/ 1.808d+01, 8.787d+00, 2.027d+00, 2.173d+00, 1.713d-01, &
                                          1.413d+01, 3.757d-01, 1.462d+01, 1.526d-01, 1.040d+01, &
                                          6.485d+00, 2.443d+00, 1.474d+01, 3.310d+01, 4.390d+02, &
                                          1.104d+02 /)
  verner_cross_sections%sig_0(16,1:16) = (/ 4.564d+04, 3.136d+02, 6.666d+00, 2.606d+00, 5.072d-04, &
                                            9.139d+00, 5.703d-01, 3.161d+01, 9.646d+03, 5.364d+01, &
                                            1.275d+01, 3.490d-01, 2.294d+04, 2.555d+01, 2.453d+01, &
                                            2.139d+02 /)
  verner_cross_sections%y_a(16,1:16) = (/ 1.000d+00, 3.442d+00, 5.454d+01, 6.641d+01, 1.986d+02, &
                                          1.656d+03, 1.460d+02, 1.611d+01, 1.438d+03, 3.641d+01, &
                                          6.583d+01, 5.411d+02, 1.529d+00, 3.821d+01, 4.405d+01, &
                                          3.288d+01 /)
  verner_cross_sections%P(16,1:16) = (/ 13.61d0, 12.81d0, 8.611d0, 8.655d0, 13.07d0, 3.626d0, &
                                        11.35d0, 8.642d0, 5.977d0, 7.09d0, 7.692d0, 7.769d0, &
                                        25.68d0, 5.037d0, 1.765d0, 2.963d0 /)
  verner_cross_sections%y_w(16,1:16) = (/ 6.385d-01, 7.354d-01, 4.109d+00, 1.863d+00, 7.880d-01, &
                                          0.000d+00, 1.503d+00, 1.153d-03, 1.492d+00, 2.310d+00, &
                                          1.678d+00, 7.033d-01, 2.738d+01, 0.000d+00, 0.000d+00, &
                                          0.000d+00 /)
  verner_cross_sections%y_0(16,1:16) = (/ 9.935d-01, 2.782d+00, 1.568d+01, 1.975d+01, 9.424d+01, &
                                          0.000d+00, 2.222d+02, 1.869d+01, 1.615d-03, 1.775d+01, &
                                          3.426d+01, 2.279d+02, 2.203d-02, 0.000d+00, 0.000d+00, &
                                          0.000d+00 /)
  verner_cross_sections%y_1(16,1:16) = (/ 0.2486d0, 0.1788d0, 9.421d0, 3.361d0, 0.6265d0, 0.d0, &
                                          4.606d0, 0.3037d0, 0.4049d0, 1.663d0, 0.137d0, 1.172d0, &
                                          0.01073d0, 0.d0, 0.d0, 0.d0 /)

  ! Iron
  verner_cross_sections%E_th(26,1:26) = (/ 7.902d+00, 1.619d+01, 3.065d+01, 5.480d+01, &
                                            7.501d+01, 9.906d+01, 1.250d+02, 1.511d+02, &
                                            2.336d+02, 2.621d+02, 2.902d+02, 3.308d+02, &
                                            3.610d+02, 3.922d+02, 4.570d+02, 4.893d+02, &
                                            1.262d+03, 1.358d+03, 1.456d+03, 1.582d+03, &
                                            1.689d+03, 1.799d+03, 1.950d+03, 2.046d+03, &
                                            8.829d+03, 9.278d+03 /)
  verner_cross_sections%E_max(26,1:26) = (/ 66.d0, 76.17d0, 87.05d0, 106.7d0, 128.8d0,  &
                                            152.7d0, 178.3d0, 205.5d0, 921.1d0, 959.d0, &
                                            998.3d0, 1039.d0, 1081.d0, 1125.d0, 1181.d0, &
                                            1216.d0, 7651.d0, 7769.d0, 7918.d0, 8041.d0, &
                                            8184.d0, 8350.d0, 8484.d0, 8638.d0, 50000.d0, &
                                            50000.d0 /)
  verner_cross_sections%E_0(26,1:26) = (/ 5.461d-02, 1.761d-01, 1.698d-01, 2.544d+01, &
                                          7.256d-01, 2.656d+00, 5.059d+00, 7.098d-02, &
                                          6.741d+00, 6.886d+01, 8.284d+00, 6.295d+00, &
                                          1.317d-01, 8.509d-01, 5.555d-02, 2.873d+01, &
                                          3.444d-01, 3.190d+01, 7.519d-04, 2.011d+01, &
                                          9.243d+00, 9.713d+00, 4.575d+01, 7.326d+01, &
                                          1.057d+03, 2.932d+02 /)
  verner_cross_sections%sig_0(26,1:26) = (/ 3.062d-01, 4.365d+03, 6.107d+00, 3.653d+02, &
                                            1.523d-03, 5.259d-01, 2.420d+04, 1.979d+01, &
                                            2.687d+01, 6.470d+01, 3.281d+00, 1.738d+00, &
                                            2.791d-03, 1.454d-01, 2.108d+02, 1.207d+01, &
                                            1.452d+00, 2.388d+00, 6.066d-05, 4.455d-01, &
                                            1.098d+01, 7.204d-02, 2.580d+04, 1.276d+01, &
                                            1.195d+01, 8.099d+01 /)
  verner_cross_sections%y_a(26,1:26) = (/ 2.671d+07, 6.298d+03, 1.555d+03, 8.913d+00, &
                                          3.736d+01, 1.450d+01, 4.850d+04, 1.745d+04, &
                                          1.807d+02, 2.062d+01, 5.360d+01, 1.130d+02, &
                                          2.487d+03, 1.239d+03, 2.045d+04, 5.150d+02, &
                                          3.960d+02, 2.186d+01, 1.606d+06, 4.236d+01, &
                                          7.637d+01, 1.853d+02, 1.358d+00, 4.914d+01, &
                                          5.769d+01, 3.288d+01 /)
  verner_cross_sections%P(26,1:26) = (/ 7.923d0, 5.204d0, 8.055d0, 6.538d0, 17.67d0, &
                                        16.32d0, 2.374d0, 6.75d0, 6.29d0, 4.111d0, &
                                        8.571d0, 8.037d0, 9.791d0, 8.066d0, 6.033d0, &
                                        3.846d0, 10.13d0, 9.589d0, 8.813d0, 9.724d0, &
                                        7.962d0, 8.843d0, 26.04d0, 4.941d0, 1.718d0, &
                                        2.963d0 /)
  verner_cross_sections%y_w(26,1:26) = (/ 2.069d+01, 1.141d+01, 8.698d+00, 5.602d-01, &
                                          5.064d+01, 1.558d+01, 2.516d-03, 2.158d+02, &
                                          2.387d-04, 2.778d-04, 3.279d-01, 3.096d-01, &
                                          6.938d-01, 4.937d-01, 1.885d-03, 0.000d+00, &
                                          1.264d+00, 2.902d-02, 4.398d+00, 2.757d+00, &
                                          1.748d+00, 9.551d-03, 2.723d+01, 0.000d+00, &
                                          0.000d+00, 0.000d+00 /)
  verner_cross_sections%y_0(26,1:26) = (/ 1.382d+02, 9.272d+01, 1.760d+02, 0.000d+00, &
                                          8.871d+01, 3.361d+01, 4.546d-01, 2.542d+03, &
                                          2.494d+01, 1.190d-05, 2.971d+01, 4.671d+01, &
                                          2.170d+03, 4.505d+02, 2.706d-04, 0.000d+00, &
                                          2.891d+01, 3.805d+01, 1.915d+06, 6.847d+01, &
                                          4.446d+01, 1.702d+02, 3.582d-02, 0.000d+00, &
                                          0.000d+00, 0.000d+00 /)
  verner_cross_sections%y_1(26,1:26) = (/ 2.481d-01, 1.075d+02, 1.847d+01, 0.000d+00, &
                                          5.280d-02, 3.743d-03, 2.683d+01, 4.672d+02, &
                                          8.251d+00, 6.570d-03, 5.220d-01, 1.425d-01, &
                                          6.852d-03, 2.504d+00, 1.628d+00, 0.000d+00, &
                                          3.404d+00, 4.805d-01, 3.140d+01, 3.989d+00, &
                                          3.512d+00, 4.263d+00, 8.712d-03, 0.000d+00, &
                                          0.000d+00, 0.000d+00 /)

END SUBROUTINE initialize_cross_sections

FUNCTION getCrosssection_rtz(lambda, element, ion) result(cross_sec)
   use constants, only: eV2erg, c_cgs, hplanck
   implicit none
   real(KIND=8), intent(in)::lambda
   integer, intent(in)::element, ion
   real(KIND=8):: cross_sec
   real(KIND=8) :: x, y, F, E

   ! Initialize cross section to 0
   cross_sec = 0.d0

   ! Convert lambda into ev
   E = hplanck * c_cgs/(lambda*1.d-8) / eV2erg         ! photon energy in ev

   ! Deal with molecular hydrogen separately
   if (element.eq.1 .and. ion.eq.3) then
      cross_sec = 0.
      if (E .gt. 11.20 .and. E .lt. 13.59) cross_sec = 2.47d-18
      if (E .ge. 13.59 .and. E .le. 15.21) cross_sec = 0.0d0
      if (E .gt. 15.21 .and. E .le. 15.45) cross_sec = 0.09d-18
      if (E .gt. 15.70 .and. E .le. 15.95) cross_sec = 1.15d-18
      if (E .gt. 15.95 .and. E .le. 16.20) cross_sec = 3.00d-18
      if (E .gt. 16.20 .and. E .le. 16.40) cross_sec = 5.00d-18
      if (E .gt. 16.40 .and. E .le. 16.65) cross_sec = 6.75d-18
      if (E .gt. 16.65 .and. E .le. 16.85) cross_sec = 8.00d-18
      if (E .gt. 16.85 .and. E .le. 17.00) cross_sec = 9.00d-18
      if (E .gt. 17.00 .and. E .le. 17.20) cross_sec = 9.50d-18
      if (E .gt. 17.20 .and. E .le. 17.65) cross_sec = 9.80d-18
      if (E .gt. 17.65 .and. E .le. 18.10) cross_sec = 10.10d-18
      if (E .gt. 18.10) cross_sec = (10.10d-18) * ((18.10d0/E)**3.0d0)
      return
   end if

   x = (E / verner_cross_sections%E_0(element,ion)) - verner_cross_sections%y_0(element,ion)
   y = sqrt( (x*x) + (verner_cross_sections%y_1(element,ion)**2.d0) )

   F = (x - 1.d0)**2
   F = F + (verner_cross_sections%y_w(element,ion)**2.d0)
   F = F * (y**(0.5d0*verner_cross_sections%P(element,ion) - 5.5d0))
   F = F * ((1.d0 + sqrt(y/verner_cross_sections%y_a(element,ion)))**(-1.d0*verner_cross_sections%P(element,ion)))

   cross_sec = verner_cross_sections%sig_0(element,ion) * F * 1.d-18
   if (E.lt.verner_cross_sections%E_th(element,ion)) then
      cross_sec = 0.d0
   endif

   if (E.gt.verner_cross_sections%E_max(element,ion)) then
      cross_sec = 0.d0
   endif

END FUNCTION getCrosssection_rtz

end module cross_sections_module