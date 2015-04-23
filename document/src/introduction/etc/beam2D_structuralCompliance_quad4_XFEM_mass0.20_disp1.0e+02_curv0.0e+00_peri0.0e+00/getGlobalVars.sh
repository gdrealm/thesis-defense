grep -r "GSOL_INTERFACE_PERIMETER:" beam2D.out | awk '{print $2;}' > GSOL_INTERFACE_PERIMETER
grep -r "GSOL_STRAIN_ENERGY:" beam2D.out | awk '{print $2;}' > GSOL_STRAIN_ENERGY
