#!/usr/bin/env python3
import sys
from datetime import datetime

# Cutoff dates for historical input changes
D_2021_07_15 = datetime(2021, 7, 15)
D_2021_08_10 = datetime(2021, 8, 10)
D_2023_01_05 = datetime(2023, 1, 5)
D_2023_01_07 = datetime(2023, 1, 7)
D_2023_01_17 = datetime(2023, 1, 17)
D_2023_02_01 = datetime(2023, 2, 1)
D_2023_02_04 = datetime(2023, 2, 4)
D_2025_10_17 = datetime(2025, 10, 17)
D_2026_05_27 = datetime(2026, 5, 27)

def parse_date(date_str):
    for fmt in ('%Y-%m-%d', '%Y/%m/%d', '%Y%m%d'):
        try:
            return datetime.strptime(date_str.strip(), fmt)
        except ValueError:
            pass
    raise ValueError(f"Could not parse date: '{date_str}'. Please use YYYY-MM-DD format.")

def get_clean_lines(filename):
    lines = []
    with open(filename, 'r') as f:
        for line in f:
            line_clean = line.strip()
            if not line_clean or line_clean.startswith(';'):
                continue
            # Remove comments starting with ';'
            if ';' in line_clean:
                line_clean = line_clean.split(';')[0].strip()
            if not line_clean:
                continue
            lines.append(line_clean)
    return lines

def parse_inp(lines, date):
    ptr = 0
    def next_line():
        nonlocal ptr
        if ptr >= len(lines):
            raise ValueError(f"Unexpected end of input file at parameter line {ptr}")
        val = lines[ptr]
        ptr += 1
        return val

    # Helper to get tokens of next line
    def get_tokens():
        line = next_line()
        return line.replace(',', ' ').split()

    # 1. Mesh parameters
    tokens = get_tokens()
    nex, nez = int(tokens[0]), int(tokens[1])
    
    tokens = get_tokens()
    x0, z0 = float(tokens[0]), float(tokens[1])
    
    tokens = get_tokens()
    rxbo, rzbo = float(tokens[0]), float(tokens[1])
    
    tokens = get_tokens()
    ircoord = int(tokens[0])
    coordfile = tokens[1] if len(tokens) > 1 else ""

    tokens = get_tokens()
    nzonx = int(tokens[0])
    zones_x = []
    for _ in range(nzonx):
        tk = get_tokens()
        zones_x.append((int(tk[0]), float(tk[1])))

    tokens = get_tokens()
    nzony = int(tokens[0])
    zones_y = []
    for _ in range(nzony):
        tk = get_tokens()
        zones_y.append((int(tk[0]), float(tk[1])))

    # 2. Tracer/Marker Section (old format only)
    if date < D_2023_01_17:
        if date < D_2021_07_15:
            # iint_marker, iint_tracer
            get_tokens()
            nzone_marker = int(get_tokens()[0])
            for _ in range(nzone_marker):
                get_tokens()
            # nzone_tracer, dt_outtracer
            get_tokens()
            for _ in range(nzone_marker):
                get_tokens()
        else:
            # iint_marker, i_junk
            get_tokens()
            nzone_marker = int(get_tokens()[0])
            for _ in range(nzone_marker):
                get_tokens()
            # n_junk, d_junk
            get_tokens()

    # 3. Mechanical boundary conditions
    tokens = get_tokens()
    nystressbc, nydrsides = int(tokens[0]), int(tokens[1])
    
    tokens = get_tokens()
    nofbc = int(tokens[0])
    bcs = []
    for _ in range(nofbc):
        tk = get_tokens()
        nofside = int(tk[0])
        nbc1 = int(tk[1])
        nbc2 = int(tk[2])
        nbc = int(tk[3])
        vals = list(map(float, tk[4:13]))
        bcs.append((nofside, nbc1, nbc2, nbc) + tuple(vals))

    tokens = get_tokens()
    nyhydro = int(tokens[0])
    pisos = float(tokens[1])
    iphsub = int(tokens[2])
    drosub = float(tokens[3])
    damp_vis = float(tokens[4])

    tokens = get_tokens()
    g = float(tokens[0])

    # 4. Thermal
    tokens = get_tokens()
    i_prestress = int(tokens[0])

    if date >= D_2025_10_17:
        tokens = get_tokens()
        extra_pres = float(tokens[0])
    else:
        extra_pres = 0.0e9

    tokens = get_tokens()
    itherm = int(tokens[0])
    
    tokens = get_tokens()
    istress_therm = int(tokens[0])
    
    tokens = get_tokens()
    ishearh = int(tokens[0])
    
    tokens = get_tokens()
    t_top = float(tokens[0])
    
    tokens = get_tokens()
    t_bot = float(tokens[0])

    tokens = get_tokens()
    hs, hr = float(tokens[0]), float(tokens[1])
    
    tokens = get_tokens()
    itemp_bc, bot_bc = int(tokens[0]), float(tokens[1])

    tokens = get_tokens()
    irtemp = int(tokens[0])
    
    tokens = get_tokens()
    tempfile = tokens[0]

    if date < D_2023_01_17:
        # time_scale
        get_tokens()

    tokens = get_tokens()
    nzone_age = int(tokens[0])
    zones_age = []
    for _ in range(nzone_age):
        first_line = get_tokens()
        if len(first_line) >= 10:  # Old 1-line format
            if date < D_2023_01_05:
                ictherm = 1
                age_1 = float(first_line[0])
                tp1 = 0.0
                tp2 = 0.0
                hcs = list(map(float, first_line[1:5]))
                iphs = list(map(int, first_line[5:10]))
                ixtb1 = int(first_line[10])
                ixtb2 = int(first_line[11])
                nph_layer = 5
            else:
                ictherm = int(first_line[0])
                age_1 = float(first_line[1])
                tp1 = float(first_line[2])
                tp2 = float(first_line[3])
                hcs = list(map(float, first_line[4:8]))
                iphs = list(map(int, first_line[8:13]))
                ixtb1 = int(first_line[13])
                ixtb2 = int(first_line[14])
                nph_layer = 5
            zones_age.append((ictherm, age_1, tp1, tp2, nph_layer, hcs, iphs, ixtb1, ixtb2))
        else:  # 3-line format
            if len(first_line) == 7:  # nph_layer on line 1
                ictherm = int(first_line[0])
                age_1 = float(first_line[1])
                tp1 = float(first_line[2])
                tp2 = float(first_line[3])
                nph_layer = int(first_line[4])
                ixtb1 = int(first_line[5])
                ixtb2 = int(first_line[6])

                line2 = get_tokens()
                hcs = list(map(float, line2))

                line3 = get_tokens()
                iphs = list(map(int, line3))
            else:  # nph_layer prepended to line 2 (newest)
                ictherm = int(first_line[0])
                age_1 = float(first_line[1])
                tp1 = float(first_line[2])
                tp2 = float(first_line[3])
                ixtb1 = int(first_line[4])
                ixtb2 = int(first_line[5])

                line2 = get_tokens()
                nph_layer = int(line2[0])
                hcs = list(map(float, line2[1:]))

                line3 = get_tokens()
                iphs = list(map(int, line3))
            
            # Pad hcs and iphs to 5 elements for uniform representation
            hcs = hcs + [0.0] * (5 - len(hcs))
            iphs = iphs + [1] * (5 - len(iphs))
            zones_age.append((ictherm, age_1, tp1, tp2, nph_layer, hcs, iphs, ixtb1, ixtb2))

    # 5. Rheology
    tokens = get_tokens()
    nphase = int(tokens[0])
    phases = []
    for _ in range(nphase):
        tk = get_tokens()
        irheol = int(tk[0])
        visc = float(tk[1])
        den = float(tk[2])
        alfa = float(tk[3])
        beta = float(tk[4])
        pln = float(tk[5])
        acoef = float(tk[6])
        eactiv = float(tk[7])
        
        # Auto-detect vactiv based on parameter count
        if len(tk) >= 25:
            vactiv = float(tk[8])
            rest = list(map(float, tk[9:25]))
        else:
            vactiv = 0.0
            rest = list(map(float, tk[8:24]))

        phases.append((irheol, visc, den, alfa, beta, pln, acoef, eactiv, vactiv) + tuple(rest))

    # 6. Phase layers / horizontal layers (only if date < 2023-01-17)
    if date < D_2023_01_17:
        mphase = int(get_tokens()[0])
        nphasl = int(get_tokens()[0])
        for _ in range(nphasl):
            get_tokens()

    # 7. Initial phase distribution
    tokens = get_tokens()
    irphase = int(tokens[0])
    phasefile = get_tokens()[0]

    # 8. Inhomogeneities
    tokens = get_tokens()
    inhom = int(tokens[0])
    inhoms = []
    for _ in range(inhom):
        tk = get_tokens()
        ix1_i = int(tk[0])
        ix2_i = int(tk[1])
        iy1_i = int(tk[2])
        iy2_i = int(tk[3])
        inphase_i = int(tk[4])
        igeom_i = int(tk[5])
        amp = float(tk[6])
        inhoms.append((ix1_i, ix2_i, iy1_i, iy2_i, inphase_i, igeom_i, amp))

    # 9. Limits & Healing
    tokens = get_tokens()
    ten_off = float(tokens[0])
    
    tokens = get_tokens()
    tau_heal = float(tokens[0])
    
    v_lims = get_tokens()
    v_min = float(v_lims[0])
    v_max = float(v_lims[1])
    ivis_shape = int(v_lims[2])
    efoldc = float(v_lims[3])

    # Melting parameters
    if date < D_2021_08_10:
        line1 = get_tokens()
        itype_melting = int(line1[0])
        nelem_serp = int(line1[1])
        prod_magma = float(line1[2])
        nelem_dike = 1
        rho_magma = 2700.0

        line2 = get_tokens()
        angle_mzone = float(line2[0])
        fmagma_max = float(line2[1])
        ratio_mantle_mzone = float(line2[3])

        line3 = get_tokens()
        latent_heat_magma = 4.2e5
        lambda_freeze = float(line3[0])
        lambda_freeze_tdep = float(line3[1])

        line4 = get_tokens()
        weaken_ratio_plastic = float(line4[0])
        weaken_ratio_viscous = float(line4[1])

    elif date < D_2023_02_01:
        line1 = get_tokens()
        itype_melting = int(line1[0])
        nelem_serp = int(line1[1])
        prod_magma = float(line1[2])
        nelem_dike = 1
        rho_magma = 2700.0

        line2 = get_tokens()
        angle_mzone = float(line2[0])
        fmagma_max = float(line2[1])
        ratio_mantle_mzone = float(line2[3])

        line3 = get_tokens()
        latent_heat_magma = 4.2e5
        lambda_freeze = float(line3[0])
        lambda_freeze_tdep = float(line3[1])

        line4 = get_tokens()
        weaken_ratio_plastic = float(line4[0])
        weaken_ratio_viscous = float(line4[1])

    elif date < D_2023_02_04:
        line1 = get_tokens()
        itype_melting = int(line1[0])
        nelem_serp = int(line1[1])
        prod_magma = float(line1[2])
        rho_magma = float(line1[3])
        nelem_dike = 1

        line2 = get_tokens()
        angle_mzone = float(line2[0])
        fmagma_max = float(line2[1])
        ratio_mantle_mzone = float(line2[3])

        line3 = get_tokens()
        latent_heat_magma = float(line3[0])
        lambda_freeze = float(line3[1])
        lambda_freeze_tdep = float(line3[2])

        line4 = get_tokens()
        weaken_ratio_plastic = float(line4[0])
        weaken_ratio_viscous = float(line4[1])

    elif date < D_2026_05_27:
        line1 = get_tokens()
        itype_melting = int(line1[0])
        nelem_serp = int(line1[1])
        prod_magma = float(line1[2])
        rho_magma = float(line1[3])
        nelem_dike = 1

        line2 = get_tokens()
        angle_mzone = float(line2[0])
        fmagma_max = float(line2[1])
        ratio_mantle_mzone = float(line2[2])

        line3 = get_tokens()
        latent_heat_magma = float(line3[0])
        lambda_freeze = float(line3[1])
        lambda_freeze_tdep = float(line3[2])

        line4 = get_tokens()
        weaken_ratio_plastic = float(line4[0])
        weaken_ratio_viscous = float(line4[1])

    else:
        line1 = get_tokens()
        itype_melting = int(line1[0])
        nelem_serp = int(line1[1])
        nelem_dike = int(line1[2])
        prod_magma = float(line1[3])
        rho_magma = float(line1[4])

        line2 = get_tokens()
        angle_mzone = float(line2[0])
        fmagma_max = float(line2[1])
        ratio_mantle_mzone = float(line2[2])

        line3 = get_tokens()
        latent_heat_magma = float(line3[0])
        lambda_freeze = float(line3[1])
        lambda_freeze_tdep = float(line3[2])

        line4 = get_tokens()
        weaken_ratio_plastic = float(line4[0])
        weaken_ratio_viscous = float(line4[1])

    # 10. Remeshing
    rem_line = get_tokens()
    ny_rem = int(rem_line[0])
    mode_rem = int(rem_line[1])
    ntest_rem = int(rem_line[2])
    angle_rem = float(rem_line[3])

    dx_rem = float(get_tokens()[0])
    
    tokens = get_tokens()
    topo_kappa, fac_kappa = float(tokens[0]), float(tokens[1])

    # 11. Process control
    idt_scale = int(get_tokens()[0])
    
    tokens = get_tokens()
    dt_scale, tolerance = float(tokens[0]), float(tokens[1])
    
    tokens = get_tokens()
    i_rey, xReyn = int(tokens[0]), float(tokens[1])

    if_rmasses = int(get_tokens()[0])
    if_imasses = int(get_tokens()[0])
    if_visc = int(get_tokens()[0])
    if_avgsr = int(get_tokens()[0])

    tokens = get_tokens()
    frac_elastic, frac_maxwell = float(tokens[0]), float(tokens[1])
    
    tokens = get_tokens()
    movegrid, ndim = int(tokens[0]), int(tokens[1])
    
    tokens = get_tokens()
    damping, mix1, mix2 = float(tokens[0]), int(tokens[1]), int(tokens[2])

    # 12. Output
    time_max = float(get_tokens()[0])
    dtout_screen = float(get_tokens()[0])
    dtout_file = float(get_tokens()[0])

    if date < D_2023_01_17:
        nprofil = int(get_tokens()[0])
        for _ in range(nprofil):
            get_tokens()

    io_flags = list(map(int, get_tokens()))
    lastout = int(get_tokens()[0])
    dtsave_file = float(get_tokens()[0])
    lastsave = int(get_tokens()[0])

    return {
        'nex': nex, 'nez': nez, 'x0': x0, 'z0': z0, 'rxbo': rxbo, 'rzbo': rzbo,
        'ircoord': ircoord, 'coordfile': coordfile, 'nzonx': nzonx, 'zones_x': zones_x,
        'nzony': nzony, 'zones_y': zones_y,
        'nystressbc': nystressbc, 'nydrsides': nydrsides, 'nofbc': nofbc, 'bcs': bcs,
        'nyhydro': nyhydro, 'pisos': pisos, 'iphsub': iphsub, 'drosub': drosub, 'damp_vis': damp_vis,
        'g': g, 'i_prestress': i_prestress, 'extra_pres': extra_pres, 'itherm': itherm,
        'istress_therm': istress_therm, 'ishearh': ishearh, 't_top': t_top, 't_bot': t_bot,
        'hs': hs, 'hr': hr, 'itemp_bc': itemp_bc, 'bot_bc': bot_bc, 'irtemp': irtemp,
        'tempfile': tempfile, 'nzone_age': nzone_age, 'zones_age': zones_age,
        'nphase': nphase, 'phases': phases, 'irphase': irphase, 'phasefile': phasefile,
        'inhom': inhom, 'inhoms': inhoms, 'ten_off': ten_off, 'tau_heal': tau_heal,
        'v_min': v_min, 'v_max': v_max, 'ivis_shape': ivis_shape, 'efoldc': efoldc,
        'itype_melting': itype_melting, 'nelem_serp': nelem_serp, 'nelem_dike': nelem_dike,
        'prod_magma': prod_magma, 'rho_magma': rho_magma, 'angle_mzone': angle_mzone,
        'fmagma_max': fmagma_max, 'ratio_mantle_mzone': ratio_mantle_mzone,
        'latent_heat_magma': latent_heat_magma, 'lambda_freeze': lambda_freeze,
        'lambda_freeze_tdep': lambda_freeze_tdep, 'weaken_ratio_plastic': weaken_ratio_plastic,
        'weaken_ratio_viscous': weaken_ratio_viscous, 'ny_rem': ny_rem, 'mode_rem': mode_rem,
        'ntest_rem': ntest_rem, 'angle_rem': angle_rem, 'dx_rem': dx_rem, 'topo_kappa': topo_kappa,
        'fac_kappa': fac_kappa, 'idt_scale': idt_scale, 'dt_scale': dt_scale, 'tolerance': tolerance,
        'i_rey': i_rey, 'xReyn': xReyn, 'ifreq_rmasses': if_rmasses, 'ifreq_imasses': if_imasses,
        'ifreq_visc': if_visc, 'ifreq_avgsr': if_avgsr, 'frac_elastic': frac_elastic,
        'frac_maxwell': frac_maxwell, 'movegrid': movegrid, 'ndim': ndim, 'damping': damping,
        'mix1': mix1, 'mix2': mix2, 'time_max': time_max, 'dtout_screen': dtout_screen,
        'dtout_file': dtout_file, 'io_flags': io_flags, 'lastout': lastout,
        'dtsave_file': dtsave_file, 'lastsave': lastsave
    }

def write_new_inp(data, filename):
    with open(filename, 'w') as f:
        f.write("; -*- LISP -*-\n;\n; Automatically upgraded to newest GeoFLAC version\n;\n")
        
        # 1. Mesh
        f.write(";==================================================================\n")
        f.write(";             M e s h    P a r a m e t e r s \n")
        f.write(";==================================================================\n")
        f.write(f"{data['nex']},{data['nez']}            number of _elements_ in X and Z directions\n")
        f.write(f"{data['x0']},{data['z0']}           x0,z0 begin.coord\n")
        f.write(f"{data['rxbo']},{data['rzbo']}    rxbo,rzbo (size of the region)\n;\n")
        f.write(f"{data['ircoord']}, {data['coordfile']}    ircoord, coordfile: read init. coordinates from file\n;\n")
        
        f.write(f"; X direction\n")
        f.write(f"{data['nzonx']}     Number zones X-direction (0 - regular grid)\n")
        for nelz, sizez in data['zones_x']:
            f.write(f" {nelz}  {sizez}\n")
        f.write(";\n")
        
        f.write(f"; Z direction\n")
        f.write(f"{data['nzony']}     Number zones Z-direction (0 - regular grid)\n")
        for nelz, sizez in data['zones_y']:
            f.write(f" {nelz}  {sizez}\n")
        f.write(";\n")
        
        # 2. Mechanical boundary conditions
        f.write(";===================================================================\n")
        f.write(";        C o n d i t i o n s:  M e c h a n i c a l\n")
        f.write(";===================================================================\n")
        f.write(f"{data['nystressbc']}  {data['nydrsides']}    nystressbc, nydrsides\n")
        f.write(f"{data['nofbc']}    nofbc: Number of boundary conditions\n")
        f.write(";nofside  nbc1 nbc2  nbc   a       b    c     d     e     f      g     h     i \n")
        for bc in data['bcs']:
            line_str = "  ".join(f"{v:g}" if isinstance(v, float) else str(v) for v in bc)
            f.write(f"{line_str}\n")
        f.write(";\n")
        
        f.write("; Hydrostatic pressure applied at bottom\n")
        f.write(f"{data['nyhydro']}  {data['pisos']}  {data['iphsub']}  {data['drosub']}  {data['damp_vis']}\n")
        f.write(";\n")
        f.write(f"; Gravity\n{data['g']}\n")
        
        # 3. Thermal
        f.write(";=============================================================\n")
        f.write(";            C o n d i t i o n s : T h e r m a l\n")
        f.write(";=============================================================\n")
        f.write(f"{data['i_prestress']}       -iprestress: allow topo build up by isostasy\n")
        f.write(f"{data['extra_pres']}    -extra pressure (GPa)\n")
        f.write(f"{data['itherm']}        -itherm (1-mech+therm, 2-no mech)\n")
        f.write(f"{data['istress_therm']}        -istress_therm: Add thermal stresses\n")
        f.write(f"{data['ishearh']}        -ishearh: Add shear heating\n")
        f.write(f"{data['t_top']}       -t_top (Surface temperature)\n")
        f.write(f"{data['t_bot']}     -t_bot (Bottom temperature)\n;\n")
        
        f.write(f"; Radiogenic heating\n")
        f.write(f"{data['hs']}, {data['hr']}      -hs (W/kg), hr (km)\n;\n")
        f.write(f"; Bottom Boundary condition flag and value\n")
        f.write(f"{data['itemp_bc']} {data['bot_bc']}\n;\n")
        
        f.write(f"; Predefined distributions\n")
        f.write(f"{data['irtemp']}              irtemp (0,1)\n")
        f.write(f"{data['tempfile']}\n")
        f.write(f"{data['nzone_age']}              - nzone_age\n")
        
        # Write nzone_age columns in newest 3-line format
        for ictherm, age_1, tp1, tp2, nph_layer, hcs, iphs, ixtb1, ixtb2 in data['zones_age']:
            f.write(f" {ictherm}, {age_1}, {tp1}, {tp2}, {ixtb1}, {ixtb2}\n")
            hc_line = " ".join(map(str, hcs[:nph_layer-1]))
            f.write(f"    {nph_layer}, {hc_line}\n")
            iph_line = " ".join(map(str, iphs[:nph_layer]))
            f.write(f"    {iph_line}\n")
        f.write(";\n")
        
        # 4. Rheology
        f.write(";===================================================================\n")
        f.write(";                     R h e o l o g y\n")
        f.write(";===================================================================\n")
        f.write(f"{data['nphase']}  Number of Different Rheologies\n;\n")
        f.write(";irheol,_,den, alfa,  beta,    n,       A,       E,    V,    Lame:rl, Lame:rm,pls1,pls2,fric1,fric2, coh1, coh2,dilat1,dilat2,cond,    cp,     Ts,     Tl,     Tk, fk\n")
        for p in data['phases']:
            line_str = ", ".join(f"{v:g}" if isinstance(v, float) else str(v) for v in p)
            f.write(f"{line_str}\n")
        f.write(";\n")
        
        # 5. Phase distribution & Inhomogeneities
        f.write(f"; INITIAL PHASE DISTRIBUTION\n")
        f.write(f"{data['irphase']}              ; irphase\n")
        f.write(f"{data['phasefile']}\n;\n")
        
        f.write(f"; Initial heterogeneities\n")
        f.write(f"{data['inhom']}  - inhom(number of heterogeneities)\n")
        for inh in data['inhoms']:
            line_str = "  ".join(f"{v:g}" if isinstance(v, float) else str(v) for v in inh)
            f.write(f"  {line_str}\n")
        f.write(";\n")
        
        # 6. Mechanical limits & Melting
        f.write(f"; Tension cut off\n{data['ten_off']}\n;\n")
        f.write(f"; linear healing parameter\n{data['tau_heal']}\n;\n")
        f.write(f"; VISCOSITY LIMIT\n")
        f.write(f"{data['v_min']}, {data['v_max']}, {data['ivis_shape']}, {data['efoldc']}\n;\n")
        
        f.write(f"; Magma:\n")
        f.write(f"; itype_melting, nelem_serp, nelem_dike, prod_magma, rho_magma\n")
        f.write(f"{data['itype_melting']}, {data['nelem_serp']}, {data['nelem_dike']}, {data['prod_magma']}, {data['rho_magma']}\n")
        f.write(f"; angle_mzone, fmagma_max, ratio_mantle_mzone\n")
        f.write(f"{data['angle_mzone']}, {data['fmagma_max']}, {data['ratio_mantle_mzone']}\n")
        f.write(f"; latent_heat_magma, lambda_freeze, lambda_freeze_tdep\n")
        f.write(f"{data['latent_heat_magma']}, {data['lambda_freeze']}, {data['lambda_freeze_tdep']}\n")
        f.write(f"; weaken_ratio_plastic, weaken_ratio_viscous\n")
        f.write(f"{data['weaken_ratio_plastic']}, {data['weaken_ratio_viscous']}\n;\n")
        
        # 7. Remeshing
        f.write(";=================================================================\n")
        f.write(";                       R e m e s h i n g\n")
        f.write(";=================================================================\n")
        f.write(f"{data['ny_rem']}  {data['mode_rem']}  {data['ntest_rem']}  {data['angle_rem']}\n;\n")
        f.write(f"; dx_rem\n{data['dx_rem']}\n;\n")
        f.write(f"; topo diffusivity and amplification factor\n{data['topo_kappa']}  {data['fac_kappa']}\n;\n")
        
        # 8. Process control
        f.write(";=================================================================\n")
        f.write(";                   P r o c e s s   c o n t r o l\n")
        f.write(";=================================================================\n")
        f.write(f"{data['idt_scale']}         - idt_scale\n")
        f.write(f"{data['dt_scale']},{data['tolerance']}  - dt_scale, tolerance\n")
        f.write(f"{data['i_rey']},{data['xReyn']}   - i_rey, Reyn\n")
        f.write(f"{data['ifreq_rmasses']}       - frequency of evaluation of real masses\n")
        f.write(f"{data['ifreq_imasses']}       - frequency of evaluation of inertial masses\n")
        f.write(f"{data['ifreq_visc']}       - frequency of evaluation of Non-Newtonian visc\n")
        f.write(f"{data['ifreq_avgsr']}       - frequency of averaging strain rate\n;\n")
        
        f.write(f"{data['frac_elastic']},{data['frac_maxwell']}\n")
        f.write(f"{data['movegrid']},{data['ndim']}\n")
        f.write(f"{data['damping']},{data['mix1']},{data['mix2']}\n;\n")
        
        # 9. Output
        f.write(";======================================================================\n")
        f.write(";                             O U T P U T\n")
        f.write(";======================================================================\n")
        f.write(f"{data['time_max']}  ; max time in kyrs\n")
        f.write(f"{data['dtout_screen']}  ; dtout_screen in kyrs\n")
        f.write(f"{data['dtout_file']}  ; dtout_file in kyrs\n;\n")
        
        f.write("; Variables to print\n")
        io_str = "   ".join(map(str, data['io_flags']))
        f.write(f"  {io_str}\n;\n")
        
        f.write(f"; output for last step only (1) or each nout step (0)\n")
        f.write(f"{data['lastout']}\n;\n")
        f.write(f"; Time interval for process saving\n")
        f.write(f"{data['dtsave_file']}\n;\n")
        f.write(f"; saving last step only (1) or each nsave step (0)\n")
        f.write(f"{data['lastsave']}\n")

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 update_inp.py <input_filename> <YYYY-MM-DD_code_date>")
        sys.exit(1)
        
    filename = sys.argv[1]
    date_str = sys.argv[2]
    
    try:
        date = parse_date(date_str)
    except ValueError as e:
        print(e)
        sys.exit(1)
        
    print(f"Tokenizing and parsing '{filename}' targeting code date {date.strftime('%Y-%m-%d')}...")
    lines = get_clean_lines(filename)
    
    try:
        data = parse_inp(lines, date)
    except Exception as e:
        print(f"Error parsing input file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
        
    # Write to a new file named with '_upgraded.inp' suffix
    out_filename = filename
    if filename.endswith('.inp'):
        out_filename = filename[:-4] + '_upgraded.inp'
    else:
        out_filename = filename + '_upgraded.inp'
        
    print(f"Writing upgraded input file to '{out_filename}'...")
    write_new_inp(data, out_filename)
    print("Upgrade completed successfully!")

if __name__ == '__main__':
    main()
