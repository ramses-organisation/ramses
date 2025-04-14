import numpy as np
import math


def check_solution(
    data,
    test_name,
    tolerance=None,
    threshold=2.0e-14,
    norm_min=1.0e-30,
    min_variance=1.0e-14,
    overwrite=False,
):
    """
    Check results against reference solution

    Parameters
    ----------
    data : dict
        Dictionary containing the data to be checked
    test_name : str
        Name of the test to be checked
    tolerance : dict, optional
        Dictionary containing the allowed relative difference between the sum over all
        cells and reference value. The default is None.
        Example: tolerance={"density":1.0e-14}.
    threshold : float, optional
        Relative value below which a vector component is set to zero.
        The default is 2.0e-14.
    norm_min : float, optional
        Minimum value for norm, to protect against null vectors. The default is 1.0e-30.
    min_variance : float, optional
        If the data differs by less than this value from the average value, it is set
        to the average. The default is 1.0e-14.
    overwrite : bool, optional
        If True, overwrite the reference solution file. The default is False.
    """

    var_tol = {"all": 3.0e-13}
    try:
        for key in tolerance.keys():
            var_tol[key] = tolerance[key]
    except AttributeError:
        pass

    # Write dummy file to avoid latex errors
    tex_file = open(test_name + ".tex", "w")
    tex_file.write(" \n")
    tex_file.close()

    # Find vectors and normalize components
    norms = dict()
    permutations = {"_x": ["_y", "_z"], "_y": ["_x", "_z"], "_z": ["_x", "_y"]}
    for key in sorted(data.keys()):
        norms[key] = 1.0
        if key.endswith("_x") or key.endswith("_y") or key.endswith("_z"):
            rawkey = key[:-2]
            suffix = key[-2:]
            ok = True
            try:
                test = len(data[rawkey + permutations[suffix][0]])
            except KeyError:
                ok = False
            try:
                test = len(data[rawkey + permutations[suffix][1]])
            except KeyError:
                ok = False
            if ok:
                norms[key] = np.sqrt(
                    data[key] ** 2
                    + data[rawkey + permutations[suffix][0]] ** 2
                    + data[rawkey + permutations[suffix][1]] ** 2
                )
                indices = norms[key] < norm_min
                norms[key][indices] = norm_min

    # Compute solution sums
    nvar = len(data.keys())
    sol = dict()
    ivar = 0
    for key in sorted(data.keys()):
        # Filter out values that are close to the average
        # This is useful if there is a constant non-zero pressure
        # with noise around in the average.
        keyAv = np.average(data[key])
        if keyAv == 0.0:
            keyData = data[key]
        else:
            keyData = np.where(
                np.abs(data[key] - keyAv) / abs(keyAv) < min_variance, keyAv, data[key]
            )
        # Perform sum
        if (
            key == "density"
            or key == "pressure"
            or key == "total_energy"
            or key == "temperature"
            or key.startswith("radiative_energy")
        ):
            solution = np.log10(np.abs(keyData))
        else:
            solution = np.where(
                np.abs(keyData) < threshold * norms[key], 0.0, np.abs(keyData)
            )

        try:
            sol[key] = math.fsum(solution)
        except TypeError:
            sol[key] = solution

    # Overwrite reference solution =====================
    if overwrite:
        print("WARNING! Over-writing reference solution")
        ref_file = open(test_name + "-ref.dat", "w")
        for key in sorted(data.keys()):
            ref_file.write("%s : %.16e\n" % (key, sol[key]))
        ref_file.close()
    # ==================================================

    # Read reference solution
    ref = dict()
    with open(test_name + "-ref.dat") as f:
        content = f.readlines()
    for line in content:
        sp = line.split(":")
        if len(sp) > 1:
            ref[sp[0].strip()] = eval(sp[1].strip())

    ok = True

    # Checking for errors
    if ref.keys() != sol.keys():
        print("The current and reference solutions do not have the same variables")
        ok = False

    # Write error table to tex file
    tex_file = open(test_name + ".tex", "w")
    # tex_file.write("\documentclass[12pt]{article}\n")
    # tex_file.write("\usepackage{graphicx,color}\n")
    # tex_file.write("\usepackage[colorlinks=true,linkcolor=blue]{hyperref}\n")
    # tex_file.write("\\begin{document}\n")
    tex_file.write("\\begin{table}[ht]\n")
    tex_file.write("\\scriptsize\n")
    tex_file.write("\\centering\n")
    tex_file.write("\\caption{" + test_name + " error summary}\n")
    tex_file.write("\\begin{tabular}{|l|l|l|l|l|}\n")
    tex_file.write("\\hline\n")
    tex_file.write("Variable & This run & Reference & Error & Tolerance\\\\\n")
    tex_file.write("\\hline\n")

    all_keys = dict()
    for key in ref.keys():
        all_keys[key] = 1
    for key in sol.keys():
        all_keys[key] = 1

    # Compute errors
    for key in sorted(all_keys.keys()):
        try:
            tol = var_tol[key]
        except KeyError:
            tol = var_tol["all"]

        try:
            this_sol = sol[key]
        except KeyError:
            this_sol = None
        try:
            this_ref = ref[key]
        except KeyError:
            this_ref = None

        if this_sol is not None and this_ref is not None:
            if this_sol == this_ref == 0.0:
                error = 0.0
            elif this_sol == 0.0 or this_ref == 0.0:
                error = np.inf
            else:
                error = abs(this_sol - this_ref) / min(abs(this_sol), abs(this_ref))
        else:
            error = np.inf

        if error > tol:
            ok = False
            output = "\\textcolor{red}{%s} & " % key.replace("_", " ")
            if this_sol is None:
                output += "\\textcolor{red}{-} & "
            else:
                output += "\\textcolor{red}{%.16e} & " % this_sol
            if this_ref is None:
                output += "\\textcolor{red}{-} & "
            else:
                output += "\\textcolor{red}{%.16e} & " % this_ref
            output += "\\textcolor{red}{%.16e} & \\textcolor{red}{%.16e} \\\\\n" % (
                error,
                tol,
            )
        else:
            output = "%s & " % key.replace("_", " ")
            if this_sol is None:
                output += "- & "
            else:
                output += "%.16e & " % this_sol
            if this_ref is None:
                output += "- & "
            else:
                output += "%.16e & " % this_ref
            output += "%.16e & %.16e\\\\\n" % (error, tol)
        tex_file.write(output)

    tex_file.write("\\hline\n")
    tex_file.write("\\end{tabular}\n")
    tex_file.write("\\end{table}\n")
    # tex_file.write("\\end{document}\n")
    tex_file.close()

    # Print message if successful
    if ok:
        print("PASSED")

    return
