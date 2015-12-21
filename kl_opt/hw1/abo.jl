#
# estimate ABO allele frequencies
#

# data
xa = 186;
xb = 38;
xab = 13;
xo = 284;

# initialization
pa = 1/3;
pb = 1/3;
po = 1/3;

# update
for i=1:10
    x = xa+xb+xab+xo;
    xnaa = xa*(pa*pa)/(pa*pa+2*pa*po);
    xnao = xa*(2*pa*po)/(pa*pa+2*pa*po);
    xnbb = xb*(pb*pb)/(pb*pb+2*pb*po);
    xnbo = xb*(2*pb*po)/(pb*pb+2*pb*po);
    pa_tmp = (2*xnaa+xnao+xab)/(2*x);
    pb_tmp = (2*xnbb+xnbo+xab)/(2*x);
    po_tmp = (xnao+xnbo+2*xo)/(2*x);
    pa = pa_tmp;
    pb = pb_tmp;
    po = po_tmp;
    println(i, " ", pa, " ", pb, " ", po);
end
