grep -r "GSOL_INTERFACE_CURVATURE:" beam2D.out | awk '{print $2;}' > GSOL_INTERFACE_CURVATURE
grep -r "GSOL_INTERFACE_PERIMETER:" beam2D.out | awk '{print $2;}' > GSOL_INTERFACE_PERIMETER
