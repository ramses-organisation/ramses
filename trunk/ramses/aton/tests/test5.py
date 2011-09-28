import os
import math

namelist = 'radiation/tests/test5.nml'
scale_d = 1.66e-27
scale_t = 3.156e15
scale_l = 9.258e22

def run_ramses():
    command = 'bin/ramses3d %s' % namelist
    result = os.system(command)
    assert result == 0

def get_output():
    cells_txt = 'tmp/cells.txt'
    command = 'utils/f90/amr2cell -inp output_00002 -out %s' % cells_txt
    result = os.system(command)
    assert result == 0

    points = []
    for line in file(cells_txt):
        points.append(map(float, line.strip().split()))
    return points

def find_nearest_point(points, x, y, z):
    d2min = None
    best_point = None
    for p in points:
        d2 = (x-p[0])**2 + (y-p[1])**2 + (z-p[2])**2
        if (d2min is None) or (d2 < d2min):
            d2min = d2
            best_point = p
    return best_point

def quantities(point):
    mH = 1.66e-24
    kB = 1.38062e-16
    scale_v = scale_l / scale_t
    scale_T2 = mH/kB * scale_v**2
    scale_pressure = scale_d * scale_v**2
    x, y, z, dx, icpu, ilevel, v1, v2, v3, v4, v5, v6, v7 = point
    rho = v1 * scale_d / mH
    xion = v6
    pressure = v5 * scale_pressure
    T2 = v5/v1 * scale_T2
    temperature = T2 / (1 + xion)
    return rho, xion, pressure, temperature

def constrain_point(points, pos,
                    rho_min, rho_max,
                    x_min, x_max,
                    p_min, p_max,
                    T_min, T_max):
    def expect_range(description, y, y_min, y_max):
        c = y_min <= y and y <= y_max
        label = c and 'pass' or 'fail'
        print '  %s: %.2e < %.2e < %.2e ? [%s]' % (
            description,
            y_min, y, y_max,
            label)
        return c

    rho, x, p, T = quantities(find_nearest_point(points, *pos))
    print 'Testing point at:', pos
    ok = True
    ok = ok and expect_range('mass density [cm^-3]', rho, rho_min, rho_max)
    ok = ok and expect_range('ionization fraction', x, x_min, x_max)
    ok = ok and expect_range('pressure [g/cm/s^2]', p, p_min, p_max)
    ok = ok and expect_range('temperature [K]', T, T_min, T_max)
    assert ok

def main():
    run_ramses()
    points = get_output()

    # The centre point at the photon source should be hot and fully ionized.
    midx = 4.0 / 8.0
    constrain_point(points,
                    (midx, midx, midx),
                    1e-4, 1e-3,   # rho
                    0.8, 1.0,     # x
                    1e-15, 1e-14, # p
                    1e4, 1e6)     # T

    # At the edge of the cube, away from the source, the gas shouldn't have
    # changed much from the initial conditions: T=100K and not ionized at all.
    constrain_point(points,
                    (0.0, 0.0, 1.0),
                    9e-4, 1.1e-3, # rho
                    0.0, 0.1,     # x
                    1e-17, 2e-17, # p
                    98, 102)      # T

if __name__ == '__main__':
    main()
