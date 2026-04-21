# CAC Documentation

<p align="right">Contributing Authors:</p>
<p align="right">CAC V1: Hao Chen</p>
<p align="right">CAC V2: Thanh Phan: tcphan@iastate.edu,tphan4@ncsu.edu
Rigelesaiyin Ji:rji@iastate.edu
Yipeng Peng: ypeng@iastate.edu</p>
<p align="right">CAC V3: Thanh Phan: tphan4@ncsu.edu
Chunxu Yan: ychunxu@ncsu.edu
</p>

<p align="right">Apr 21, 2026</p>

## Table of Contents

-   [Introduction](#introduction)
-   [Important Notes](#important-notes)
-   [Commands in Stable Version](#commands-in-stable-version)
    -   [add_atoms](#add_atoms)
    -   [add_etype](#add_etype)
    -   [atom_style](#atom_style)
    -   [atomic_strain](#atomic_strain)
    -   [balance](#balance)
    -   [boundary](#boundary)
    -   [box](#box)
    -   [change_box](#change_box)
    -   [charge](#charge)
    -   [clear](#clear)
    -   [create_atoms](#create_atoms)
    -   [create_box](#create_box)
    -   [create_elements](#create_elements)
    -   [coarse_graining](#coarse_graining)
    -   [compute centroid/stress/atom](#compute-centroidstressatom)
    -   [compute cna/atom](#compute-cnaatom)
    -   [compute displace/atom](#compute-displaceatom)
    -   [compute ids/atom](#compute-idsatom)
    -   [compute ke](#compute-ke)
    -   [compute msd](#compute-msd)
    -   [compute pe](#compute-pe)
    -   [compute pressure](#compute-pressure)
    -   [compute stress/atom](#compute-stressatom)
    -   [compute stress/mech](#compute-stressmech)
    -   [compute temp](#compute-temp)
    -   [compute vacf](#compute-vacf)
    -   [delete_atoms](#delete_atoms)
    -   [delete_elements](#delete_elements)
    -   [dimension](#dimension)
    -   [dielectric](#dielectric)
    -   [disc_elements](#disc_elements)
    -   [displace](#displace)
    -   [dump](#dump)
    -   [echo](#echo)
    -   [element_modify](#element_modify)
    -   [element_style](#element_style)
    -   [fix adaptive](#fix-adaptive)
    -   [fix addforce](#fix-addforce)
    -   [fix aveforce](#fix-aveforce)
    -   [fix ave/time](#fix-avetime)
    -   [fix balance](#fix-balance)
    -   [fix deform](#fix-deform)
    -   [fix enforce2d](#fix-enforce2d)
    -   [fix langervin](#fix-langervin)
    -   [fix momentum](#fix-momentum)
    -   [fix move](#fix-move)
    -   [fix nve](#fix-nve)
    -   [fix print](#fix-print)
    -   [fix press/berendsen](#fix-pressberendsen)
    -   [fix setforce](#fix-setforce)
    -   [fix temp/press](#fix-temppress)
    -   [fix temp/rescale](#fix-temprescale)
    -   [fix viscous](#fix-viscous)
    -   [group](#group)
    -   [jump](#jump)
    -   [label](#label)
    -   [lattice](#lattice)
    -   [log](#log)
    -   [mass](#mass)
    -   [modify_elements](#modify_elements)
    -   [newton](#newton)
    -   [neighbor](#neighbor)
    -   [neigh_modify](#neigh_modify)
    -   [next](#next)
    -   [pair_style](#pair_style)
    -   [pair_coeff](#pair_coeff)
    -   [processors](#processors)
    -   [read_data](#read_data)
    -   [region](#region)
    -   [reset_ids](#reset_ids)
    -   [reset_timestep](#reset_timestep)
    -   [run](#run)
    -   [set_atoms](#set_atoms)
    -   [set_elements](#set_elements)
    -   [thermo](#thermo)
    -   [thermo_style](#thermo_style)
    -   [timestep](#timestep)
    -   [quit](#quit)
    -   [units](#units)
    -   [variable](#variable)
    -   [velocity](#velocity)
    -   [write_data](#write_data)
    -   [write_tecplot](#write_tecplot)

## Introduction

CAC software package is build based on LAMMPS and hence the structure is similar
 to LAMMPS with some minor differences (See LAMMPS developer guide for more
   details). If something is similar to LAMMPS, see LAMMPS docmentation for
   syntax.

Current stable version: 3.1.2 21Apr2026

## Important Notes

-   When atoms/elements are created/deleted, use migrate option only for the last command
    that create/delete atoms/elements. Previous commands should turn off migrate option.

## Installation

-   In the src folder, run:
    make -j8
-   -j8 option means using 8 processors to compile in parallel for faster compilation

## Optional Packages

-   To install optional packages, run:
    make yes-package-name
-   To uninstall optional packages, run:
    make no-package-name
-   List of Optional Packages:
    - TECPLOT
    - DEEPMD

## Utilities

-   Utilities code are provided in utils folder. Check README.md file in each folder for details.

## Commands in stable version

### add_atoms

-   adding extra atom layers to a side of elements (often to create a dislocation)
    new atoms will be placed on cell positions in elements.
-   command syntax:

        add_atoms group-ID dim nlayers
        -   group-ID = ID of group of elements to add extra layers of atoms
        -   dim = +x -x +y -y +z - z, element face to add extra layers of atoms
          (in element natural coordinate system)
        -   nlayers = number of layers to add.

<a href="#table-of-contents">Back to top</a>

### add_etype

-   define an element type (must be > netypes)
-   command syntax:

        add_etype type_id element_shape apc ulx uly ulz intx inty intz
        -   type_id = index of element type
        -   element_shape = shape of elements 
            (available shapes: Quad, Tri, Hex, Tet, Oct, Pyr, Wedge)
        -   apc = number of "atoms per unit cells" 
        -   ulx,uly,ulz = size of elements (number of unit cells) 
            along x, y, z directions
        -   intx,inty,intz = number of integration points 
            (for Hex elements, set all to 4)

<a href="#table-of-contents">Back to top</a>

### atom_style

-   Define atom style
-   currently only atomic style is the default style
-   other atom styles can be added to the code in the future.

<a href="#table-of-contents">Back to top</a>

### atomic_strain
-   Perform strain calculation (same as OVITO's atomic strain modifier)
-   command syntax:

        atomic_strain  input_file  output_file keywords values
        zero or more keyword/values pairs maybe appended
        keywords = out/format or cutoff or out/config or def/gradient or element/info
        or rotation or nonaff/squared/disp or stretch/tensor or strain/tensor or wrap
        or average or ref/file or eweight or no/atom or no/element or region or step
        -   out/format = plt/szplt/ascii/atom (default = plt)
        -   cutoff = cutoff distance for strain calculation (default = 6)
        -   out/config = current/refenrece, output the results in current or reference
            configuration (default = current)
        -   def/gradient = yes/no, output 9 components of the deformation gradient 
            tensor (default = no)
        -   element/info = yes/no, output element information (default = no)
        -   rotation = yes/no, output rotation as 4 quaternions (default = no)
        -   nonaff/squared/disp = yes/no, output nonaffine squared displacement 
            (default = no)
        -   stretch/tensor = yes/no, output 6 component of the stretch tensor   
            (default = no)
        -   strain/tensor = yes/no, output 6 component of the stretch tensor 
            (default = no)
        -   wrap = yes/no, wrap positions of atoms/elements across periodic boundaries
            (default = yes)
        -   average = yes/no, average nodal position (default = yes)
        -   eweight = yes/no, element weight for rcb rebalancing workload
        -   no/atom = yes/no, exclude atoms in output
        -   no/element = yes/no, exclude elements in output
        -   region = region-ID, only output atoms/elements within the region
        -   step = timestep of this frame (used for indexing the timestep in Zone Title for Tecplot files)
        Von Mises strain and Volumetric strain are output by default

<a href="#table-of-contents">Back to top</a>

### balance

-   Similar to LAMMPS
-   atoms have uniform weights of 1.
-   element weight is number of integration points by default, can be modified
    by setting eweight option

<a href="#table-of-contents">Back to top</a>

### boundary

-   define boundary in x y z direction, similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### box

-   Similar to LAMMPS.

<a href="#table-of-contents">Back to top</a>

### change_box

-   Similar to LAMMPS.

<a href="#table-of-contents">Back to top</a>

### charge

-   define charge for each atomic type (chemical type)
-   command syntax:
    
        charge    atom_type  charge_value
        - atom_type = atomic type to assign charge value
        - charge_value = self-explanatory 

-   examples:
        charge 1 1.84 #Sr
        charge 2 2.36 #Ti
        charge 3 -1.40 #O
    
<a href="#table-of-contents">Back to top</a>

### clear

-   Clear everything in CAC

<a href="#table-of-contents">Back to top</a>

### create_atoms

-   Similar to LAMMPS.

<a href="#table-of-contents">Back to top</a>

### create_box

-   Similar to LAMMPS, no option keywords.

<a href="#table-of-contents">Back to top</a>

### create_elements

-   Command syntax:

        create_elements etype ctype style args keyword values
        -   ctype = chemical type (atom type)
        -   etype = element type (must be defined by add_etype command beforehand)
        -   style = box or region or single
            -   box args = none
            -   region args = region-ID
                region-ID = elements will only be created if center contained in the
                region
            -   single args = x y z
                x y z = coordinate of center of a single elements, node coords will be
                created according to lattice command (a lattice must be defined)
        zero or more keyword/values pairs maybe appended
        keywords = similar to LAMMPS create_atoms command (no mol keyword)

<a href="#table-of-contents">Back to top</a>

### coarse_graining

-   Coarse graining atomisitic systems into elements
-   Command syntax:

        coarse_graining  group_ID  cutoff  structure_type  element_etype  atomic_type  keywords values
        -   group_ID = group selected for coarse graining
        -   cutoff = cutoff to find first nearest neighbor
        -   structure_type = fcc or bcc or cubic/diamond
            currently, it only works for fcc structure
        -   element_type = element type of newly created elements, should be pre-defined
        -   atomic_type = atomic type of newly created elements
        zero or more keyword/values pairs maybe appended
        keywords = compress or migrate or wrap or tol or a1 or a2 or a3 or centroid
        -   compress = yes/no (default = yes)
        -   migrate = yes/no (default = yes)
        -   wrap = yes/no/x/y/z/xy/xz/yz, specify which directions to perform wrapping
            (default = xyz)
        -   tol = tolerance for aligning unit cells (default = 1e-3)
        -   a1 = a1 vector of unit cell (default from lattice command)
        -   a2 = a2 vector of unit cell (default from lattice command)
        -   a3 = a3 vector of unit cell (default from lattice command)
        -   centroid = filename, read centroid position for each grain
        -   seeds = N id_1 id_2 ... id_N, the IDs of seed atoms to start coarse graining process (ideally one for each grain)

-   example command:

        coarse_graining  group all  Cu/fcc  fcc  1  1  tol 0.01 wrap yes &
                                    a1 0 0 1 &
                                    a2 0.8660254 0 0.5 &
                                    a3 0.2886751 0.81649658 0.5

<a href="#table-of-contents">Back to top</a>

### compute centroid/stress/atom

-   Similar to LAMMPS.
-   thermo terms are not included yet.

<a href="#table-of-contents">Back to top</a>

### compute cna/atom

-   Similar to LAMMPS.
-   compute CNA for all atoms and "virtual atoms" (only on the surfaces) in group
-   For cutoff use "adaptive" only, traditional cutoff is no longer supported
-   structure id list:

        unknown: 0
        fcc: 1
        hcp: 2
        bcc: 3
        ico: 4
        other: 5

<a href="#table-of-contents">Back to top</a>

### compute displace/atom

-   Similar to LAMMPS.
-   no refresh option.
-   compute displacement of atoms and nodes from reference configurations (since
    this compute was issued)

<a href="#table-of-contents">Back to top</a>

### compute ids/atom

-   determine diamond structure for all atoms and interpolated atoms in group
-   Command syntax:

        compute id group_id ids/atom

-   structure id list:

        unknown: 0
        other: 1
        Cubic diamond: 2
        Hexagonal diamond: 3

<a href="#table-of-contents">Back to top</a>

### compute ke

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### compute msd 

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### compute pe

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### compute pressure

-   Similar to LAMMPS.
-   thermo terms are not included yet.

<a href="#table-of-contents">Back to top</a>

### compute stress/atom

-   Similar to LAMMPS.
-   thermo terms are not included yet.

<a href="#table-of-contents">Back to top</a>

### compute stress/mech

-   calculate mechanical stress on a mesh (only rectangular mesh support currently)
-   Command syntax:

        compute id group_id stress/mech dir nx ny nz xlo xhi ylo yhi zlo zhi keyword
        values
        -   dir = surface normal plane, can be either x, y, or z
        -   nx, ny, nz = number of "boxes" in x, y, z directions
        -   xlo, xhi, ylo, yhi, zlo, zhi = bounds of box to compute stress;
        keyword/value pairs:
        -   style value = standard (one style currently)
        -   one value = size of surface neighborlist
        -   position value = surface position ratio in box, from 0 to 1 (0.5 by
            default)
-   To output results, use dump stress/mech

<a href="#table-of-contents">Back to top</a>

### compute temp

-   Similar to LAMMPS.

<a href="#table-of-contents">Back to top</a>

### compute vacf

-   Similar to LAMMPS.

<a href="#table-of-contents">Back to top</a>

### delete_atoms

-   Similar to LAMMPS, only delete atoms, elements are not deleted in this command.
-   style = group or region

<a href="#table-of-contents">Back to top</a>

### delete_elements

-   Similar to delete_atoms, only delete elements, atoms are not deleted in this command.
-   style = group or region

<a href="#table-of-contents">Back to top</a>

### dielectric

-   similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### dimension

-   define simulation dimension

<a href="#table-of-contents">Back to top</a>

### disc_elements

-   command syntax:

        disc_elements   region region_ID

-   options:
    -   migrate = yes or no (default = yes)
    -   compress = yes or no (default = yes)
    -   style = atom or element (default = atom)
        -   atom args = none. Discretize elements to fully discrete atoms
        -   element args = original_shape final_shape extra_args
            -   original_shape = Hex or Wedge or Pyramid. 
                Current shape of element to split 
            -   final_shape = Final shape of element to split (args depend on original shape)
                Hex --> Hex or HexQuad or Wedge or PyrTet
                Wedge --> Wedge or PyrTet or Tet
                Pyramid --> Tet
            -   Too lazy to details everything here. See code for more details for the options or ask Thanh (^.^)
                
-   discretize all elements in region_ID into fully discrete atoms or splits to smaller elements with similar/different shapes
-   elements in region are deleted and new atoms/elements are created in place

<a href="#table-of-contents">Back to top</a>

### displace

-   Similar to displace_atoms command in LAMMPS

    displace group-ID style args keyword value ...
    -   group-ID = ID of group of atoms/elements to displace
    -   style = move or ramp or rotate or dislocation or multidisl or disp/file
        -   move: similar to LAMMPS
        -   ramp: similar to LAMMPS
        -   rotate: similar to LAMMPS
        -   dislocations: args = N ddim dim1 dim2 sdim b angle v theta
            -   N = number of dislocations, the next set of arguments must be provided for each dislocation
            -   ddim = dislocation line direction (x or y or z)
            -   dim1, dim2 = position of dislocation line in (yz or zx or xy)
            -   sdim = slip direction (+x/-x or +y/-y or +z/-z)
            -   b = length of Burgers vector (can be negative)
            -   angle = angle (in degrees) from positive dislocation line direction toward
                Burgers vector (positive toward positive slip direction)
            -   v = poisson ratio
            -   theta = rotation angle (in degrees) about the dislocation line axis
        -   multidisl: args = N ddim dim1 dim2 sdim offset r b angle v theta
            -   ddim,dim1,dim2,sdim,b,angle,v: same as above
            -   N = number of dislocations
            -   offset = distance between 2 dislocation (in distance units)
            -   r = offset ratio
        -   disl/loop: argsh = shape cx cy cz nx ny nz R bx by bz nu
            -   shape = The geometric shape of the loop (circle or square)
            -   cx, cy, cz = Coordinates of the center of the dislocation loop (distance units)
            -   nx, ny, nz = Normal vector to the loop's habit plane (unitless)
            -   R = Geometric scaling factor (distance units, must be positive). See shape details below.
            -   bx, by, bz = Burgers vector (distance units)
            -   nu = Poisson's ratio of the material (unitless)
        -   disp/file: args = file_name
            -   file_name = name of the external displacement file
    -   keyword pairs may be appended:
        -   keyword = units
            -   value = box or lattice
        -   keyword = group 
            -   value = atom or node or element
        -   keyword = sequence
            -   value = yes or no
        -   keyword = index/offset (for disp/file style only)
            -   value = offset values for node index (so that the actual index starts from 0)
        -   keyword = header/offset (for disp/file style only)
            -   value = number of header lines to skip

-   dislocations or multidisl style displace atoms/nodes positions according to
    theoretical
        dislocations displacement field to create full dislocation(s). dim1 and dim2
        must be placed correctly, if ddim is x, dim1, dim2 will be y and z coordinates
        of dislocation line respectively, if ddim is y, dim1, dim2 will be z and x
        coordinates of dislication line respectively. For multidisl style, the
        offset between dislocation increases/decreases exponentially with ratio r
        (dist_i = offset\*r^i)
-   disl/loop style: This command applies a localized displacement field to all atoms and continuum elements
    within the specified group to embed a discrete dislocation loop. The field accurately
    represents the Volterra displacement of a finite dislocation loop using the principles of
    isotropic linear elasticity.

    The displacement of any atom or node is calculated by superposing the displacement fields
    of the finite segments that construct the shape. The calculation relies on two primary
    mathematical components:
    1.  Finite Segment Elasticity: The elastic strain field radiating from each segment is
        calculated using Barnett's formula for triangular dislocation loops.
    2.  Solid Angle Cut: The discrete slip (the lattice cut) across the habit plane is generated by
        calculating the solid angle subtended by the loop at the observation point, utilizing the
        geometric algorithm proposed by Van Oosterom and Strackee.
    Interstitial vs. Vacancy Loops:
        The mathematical formulation assumes an interstitial loop by default. To simulate a vacancy
        loop (the removal of a disk of atoms), negate the components of the Burgers vector (bx, by, bz). 
    Loop Shape: 
    1.  Circle: The parameter R defines the exact radius of the circular loop. The loop is
        generated dynamically by discretizing the circle into 60 straight, finite segments lying
        on the habit plane.
    2.  Square: The parameter R acts as the half-width of the square, meaning the total side
        length of the square is 2R. The square's edges are generated using an internally
        determined set of orthogonal vectors lying on the habit plane.;

-   disp/file style displace atoms/nodes positions according to the values provided, the format for each line is as follow:

        ID index dx dy dz
        
        -   ID = ID of atom or element
        -   index = 0 if it is an atom, otherwise it is the index of the node within the element (nodeindex * basisindex + 1)
        -   dx, dy, dz = displacement values
                 
        The index value can be offset using the index/offset optional keyword. The file may contain headers and can be skipped using header/offset keyword.

-   example:

        displace all dislocations 2 x 0 0 +y 2.55 90 0.33 0 &
                                    x 0 5 +y 2.55 90 0.33 70.5288
        displace all multidisl 5 x 0 0 +y 100 1.1 2.55 90 0.33 0

        displace all disl/loop circle 50.0 50.0 50.0 0.0 1.0 0.0 25.0 2.55 0.0 0.0 0.33

-   NOTE: group-ID is ignored for dislocation, multidisl, and disl/loop styles

<a href="#table-of-contents">Back to top</a>

### dump

-   Similar to LAMMPS
-   style = tecplot or atom or cna/atom or ids/atom or stress/mech
    -   tecplot: output .dat/.plt file for Tecplot (require TECPLOT package installed if writing .plt files). Simulation box is visualized as an FE element in a separate zone.
        option keywords: stress or force or velocity or pe or displace or format or average
        -   stress  : args = yes/no/centroid, (default = no) add 6 components of virial stress values to output files, 
                centroid for 9 components from compute centroid/stress/atom (for DeepMD potential)
        -   force   : args = yes/no, (default = no) add 3 components of force to output files
        -   velocity: args = yes/no, (default = no) add 3 components of velocity to output files
        -   displace: args = yes/no, (default = no) add displacement values Ux,Uy,Uz,Utotal to output files
        -   pe      : args = yes/no, (default = no) add atomic potential energy to output files
        -   average : args = yes/no, (default = yes) average nodal positions/values
        -   format  : args = ascii/plt/szplt, (default = plt) define format of dump file. This option keyword is overwritten by the file name suffix (ascii = .dat, plt = .plt, szplt = .szplt)
    -   atom: output all atoms and interpolated atoms into dump file similar to LAMMPS
            format (no element)
        option keywords: stress or force or velocity or displace or pe or style or wrap
        -   stress  : args = yes/no, (default = no) add 6 components of virial stress values to output files
        -   force   : args = yes/no, (default = no) add 3 components of force to output files
        -   velocity: args = yes/no, (default = no) add 3 components of velocity to output files
        -   displace: args = yes/no, (default = no) add displacement values Ux,Uy,Uz,Utotal to output files
        -   pe      : args = yes/no, (default = no) add atomic potential energy to output files
        -   style   : args = interpolate/integration/node/center, (default = interpolate) output values at interpolated/integration/nodal/center points within elements
        -   wrap    : args = yes/no, (default = no) wrap atoms for each periodic boundary (parts of an element may be out of the box and need to be remapped)
    -   structure/atom: args = structure compute id (must be either cna or ids), followed by list of structure types to keep
        -   cna/atom: args = fcc and/or bcc and/or hcp and/or ico and/or other
        -   ids/atom: args = cubic and/or hex and/or other 
        output atoms and interpolated atoms with specified structure type
        require compute cna/atom or ids/atom command defined
    -   stress/mech: output mechanical stress values
        require compute stress/mech command
        need to specify the compute ID for stress/mech compute at the end of the
            command
    -   cac/data: output CAC data file for restarting simulation (velocities not included yet) or for strain calculation
        option keywords: reference or precision
        -   reference : args = current or data/file or off, (default = off) output both current and reference configuration for strain calculation, current means using the "current" configuration (at the moment dump command is defined) as the reference, data/file means using an external data file as the reference (the data file name is appended in the command)
        -   precision : args = high/normal, (default = normal) 

<a href="#table-of-contents">Back to top</a>

### echo

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### element_modify

-   Modify some parameters in element
-   Syntax:

        element_modify keyword values ...
        keyword =  id or map or max/ucell or max/gcell or max/sub/ucell or max/sub/elem or max/elem/change or max/surface/ucell 
                   or max/outer/ucell or max/edge/ucell or mass/style or rigid
        -   id value = yes or no
        -   map value = array or hash or yes
        -   max/ucell value = max number of interpolation atoms in each element 
            (= 4000 by default)
        -   max/gcell value = max number of integration points in each element 
            (= 64 by default)
        -   max/sub/ucell value = max number of interpolation atoms in each sub-element
            (= 500 by default)
        -   max/sub/elem value = max number of sub-element in an element
            (= 125 by default)
        -   max/elem/change value = throw error if element size is increased by this value
            (= 2.0 by default)
        -   max/surface/ucell value = max number of interpolation atoms in each element surface
            (= 400 by default)
        -   max/outer/ucell value = max number of interpolation atoms in each element outer surfaces
            (= 2400 by default)
        -   max/edge/ucell value = max number of interpolation atoms in each element edge
            (= 20 by default)
        -   nodal/force/style value = lumped or consistent
            (default = consistent)
        -   rigid value = yes or no

<a href="#table-of-contents">Back to top</a>

### element_style

-   Define element style
-   Syntax:

        element_style cac max_npe max_apc
        -   max_npe = max number of nodes per element
        -   max_apc = max number of atoms per unit cell
        option keyword: wscale
        -   wscale: args = nodew edgew surfacew innerw (specify weight ratios for node:edge:surface:inner integration points)

<a href="#table-of-contents">Back to top</a>

### fix adaptive

-   Syntax:

    fix ID group_ID adaptive N d splipplane region_ID compute_ID perfect_crystal_ID
    -   N = Nevery (must be > 0)
    -   d = width between splip planes
    -   splipplane = xy or yz or xz
    -   region_ID = ID of region that involve (only use when fix is initialized)
    -   compute_ID = ID of compute that determine crystal structure (cna/atom or ids/atom)
    -   perfect_crystal = structure ID of perfect crystal in compute

-   This fix choose out atoms and elements in group with in region at fix initialization.
-   The atoms in group and region are used to check for dislocation. Non-perfect 
        crystal atoms in 2 adjacent atomic planes are associated 
        with a dislocation in that slip plane. The elements in 
        group and region will be splitted if needed.
-   If a resulting element from splitting is a 2D element 
        (ncell = 1 in a direction) it will be discretized 
        into planar (Quad) element. The number of integration 
        points of a new element will stay the same or be at 
        most size of ncells in the splitting direction.
-   Data of old elements will be passed on to new atoms/elements through interpolation.

<a href="#table-of-contents">Back to top</a>

### fix addforce

-   Similar to LAMMPS
-   option keyword:
    -   group args = element or node
        -   element: choose nodes by elements in group
        -   node: choose nodes by nodes in group
-   default groups: all, atom, element

<a href="#table-of-contents">Back to top</a>

### fix aveforce

-   Similar to LAMMPS
-   option keyword:
    -   group args = element or node
        -   element: choose nodes by elements in group
        -   node: choose nodes by nodes in group
-   default groups: all, atom, element

<a href="#table-of-contents">Back to top</a>

### fix ave/time

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### fix balance

-   Similar to LAMMPS
-   Option keywords:
    -   eweight = weight of element (weight of one atom is 1)
        = number of integration points by default

<a href="#table-of-contents">Back to top</a>

### fix deform

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### fix enforce2d

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### fix halt

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### fix langervin

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### fix momentum

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### fix move

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### fix nve

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### fix print

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### fix press/berendsen

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### fix setforce

-   Similar to LAMMPS
-   option keyword:
    -   group args = element or node
        -   element: choose nodes by elements in group
        -   node: choose nodes by nodes in group

<a href="#table-of-contents">Back to top</a>

### fix temp/press

-   

<a href="#table-of-contents">Back to top</a>

### fix temp/rescale

-   

<a href="#table-of-contents">Back to top</a>

### fix viscous

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### group

-   Similar to LAMMPS
-   style = region or subtract or union or intersect or atom or element
-   for atom/element style
-   option keyword (for region style only):
    -   style args = center or allnode or onenode or oneatom
        determine how elements are considered to be in region
        -   center: if element center is in region (default)
        -   allnode: if all element nodes are in region
        -   onenode: if at least one element node is in region
        -   oneatom: if at least one interpolated atom is in region

<a href="#table-of-contents">Back to top</a>

### jump

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### label

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### lattice

-   Similar to LAMMPS
-   added built-in styles: pri/fcc, pri/bcc, pri/hcp, pri/cubic/diamond

<a href="#table-of-contents">Back to top</a>

### log

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### mass

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### newton

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### neighbor

-   Similar to LAMMPS, only bin style for now

<a href="#table-of-contents">Back to top</a>

### neigh_modify

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### next

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### pair_style

-   Similar to LAMMPS
-   Supported potentials:
    -   buck
    -   coul_wolf
    -   deepmd (require DEEPMD package installed)
    -   eam
    -   eam/alloy
    -   eam/fs
    -   gauss_cut
    -   hybrid
    -   hybrid/overlay
    -   lj/cut
    -   lj/mdf
    -   meam
    -   morse
    -   sw
    -   table
-   other potential style can be added by modifying the existing code from LAMMPS

<a href="#table-of-contents">Back to top</a>

### pair_coeff

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### processors

-   Similar to LAMMPS
-   style = file
    -   file arg = outfile (name of file to write 3D grid of processors to

<a href="#table-of-contents">Back to top</a>

### read_data

-   Similar to LAMMPS
-   keyword = add or offset or shift or extra/atom/types or extra/element/types or group or reference
    -   add arg = append or merge
    -   offset args = toff aoff
        -   toff = offset to add to atom types (and element ctypes)
        -   eoff = offset to add to element types
    -   shift args = Sx Sy Sz
        -   Sx,Sy,Sz = distance to shift atoms/elements when adding to system
    -   extra/atom/types arg = # of extra atom types
    -   extra/element/types arg = # of extra element types
    -   group arg = groupID
        -   groupID = add atoms/elements in data file to this group (create a new group if not exist already)
    -   reference arg = cur/frame or ref/frame or off (default = off), use when reading data file that include a reference configuration (from dump cac/data command), selecting current or reference configuration to read

<a href="#table-of-contents">Back to top</a>

### region

-   Similar to LAMMPS
-   Units are set to box by default
-   style = block or sphere or cylinder or plane or prism

<a href="#table-of-contents">Back to top</a>

### reset_ids

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### run

-   Similar to LAMMPS
-   no keyword

<a href="#table-of-contents">Back to top</a>

### set_atoms

-   Similar to LAMMPS, operating only on atoms

<a href="#table-of-contents">Back to top</a>

### set_elements

-   Similar to set_atoms command with type --> etype or ctype
-   Operating only on elements
-   New style: group/ctype, selecting elements in group and basis atoms with specific ctype to assign new ctype. Note: can only set ctype for this style
    Example: the command below will select the elements in group noforce, changing basis atoms with ctype 1 to ctype 10

    set_elements group/ctype noforce 1 ctype 10
     

<a href="#table-of-contents">Back to top</a>

### thermo

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### thermo_style

-   Similar to LAMMPS (see LAMMPS doc for keyword meanings)
-   style = one or custom
    -   one args = none
    -   custom args = list of keywords
        possible keywords: 

            step, elapsed, elaplong, dt, time, cpu, hcpu
            tpcpu, spcpu, cpuremain, hcpuremain, part, timeremain,
            imbalance, elements, atoms, press, pe, vol, lx, ly, lz,
            xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz, xlat, ylat, zlat,
            pxx, pyy, pzz, pxy, pxz, pyz, nbuild, ndanger.
            c_ID, c_ID[I], c_ID[I][j], f_ID, f_ID[I], f_ID[I][j],
            v_name

-   style one print out step cpu and press (set by default, volume colume will be added if box is changing)

<a href="#table-of-contents">Back to top</a>

### timestep

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### quit

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### units

-   define unit systems similar to LAMMPS
-   metal style now use Bar for pressure to be consistent with LAMMPS.
-   metalgpa style use GPa instead of Bar for pressure unit,

<a href="#table-of-contents">Back to top</a>

### variable

-   Similar to LAMMPS

<a href="#table-of-contents">Back to top</a>

### velocity

-   Similar to LAMMPS
-   units are set to box by default
-   option keyword:
    -   group args = element or node
        -   element: choose nodes by elements in group
        -   node: choose nodes by nodes in group

<a href="#table-of-contents">Back to top</a>

### write_data

-   Write a data file that can be read into another CAC simulation
-   if no element, the data file can be read into LAMMPS, otherwise, 
        it can't be read by LAMMPS.
-   option keywords: noinit or novelocity or noelement or noatom
    -   noinit    : yes/no (default = no)
    -   novelocity: yes/no (default = no)
    -   noelement : yes/no (default = no)
    -   noatom    : yes/no (default = no)

<a href="#table-of-contents">Back to top</a>

### write_tecplot

-   Write visualization file for TecPlot (require TECLOT package installed)
-   Supported file type: .dat, .plt, .szplt
-   option keywords: format or noinit
    -   format: ascii/plt/szplt (default = plt)
    -   noinit: yes/no (default = no)

<a href="#table-of-contents">Back to top</a>
