# Simple Laplace transform calculator
# made by 00glitched
# pre-run:
# symbolic install command:
pkg install -forge symbolic;
pkg load symbolic;
syms x;syms t;syms s;
syms n;syms m;
cont=1;
while(cont==1)
    disp("\n\nOptions:")
    disp("t : f(t) -> F(s) # Laplace transform");
    disp("s : F(s) -> f(t) # inverse Laplace");
    disp("c : f(t) *  g(t) # convolution");
    disp("h :              # examples");
    disp("e :              # exit");
    opt = input("\nSelect: ", "s");
    if(opt == "t")
        # u[n]  ->  heaviside(t-n)  and   f(t-n)
        # d[n]  ->  dirac(t-n)
        # L(ft)
        disp ("\n\n\L { f(t) } :\n");
        ft = sym (input ("\nf(t) = ", "s"));
        disp ("\nf(t) =\n"); pretty (ft, "ascii");
        #Fs = int(fts*exp(-s*t), t, 0, limit(1/t, t, 0, 'right')) #integral method
        FS  = laplace (ft);
        disp ("\nF(s) =\n"); pretty (FS, "ascii");
        disp ("\n\nLaTex code (to use in docs):\n");
        disp ("f(t):");disp(latex(ft));disp("");
        disp ("F(s):");disp(latex(FS));disp("");
    end
    if (opt == "s")
        # L⁻1(Fs)
        disp ("\n\n\L^-1 { F(s) } :\n");
        FS = sym (input ("\nF(s) = ", "s"));
        disp ("\nF(s) =\n"); pretty (FS, "ascii");
        ft  = ilaplace(FS);
        disp ("\nf(t) =\n"); pretty (ft, "ascii");
        disp ("\n\nLaTex code (to use in docs):\n");
        disp ("f(t):");disp(latex(ft));disp("");
        disp ("F(s):");disp(latex(FS));disp("");
    end
    if (opt == "c")
        # convolution
        disp ("\n\nf * g :\n");
        f = sym (input ("\nf(x) = ", "s"));
        disp ("\nf(x) =\n"); pretty (f, "ascii");
        fxt= subs(f,x,x-t);
        disp ("\nf(x-t) =\n"); pretty (fxt, "ascii");
        g  = sym (input ("\ng(x) = ", "s"));
        disp ("\ng(x) =\n"); pretty (g, "ascii");
        fg = int (fxt*g, x, 0, t);
        disp ("\n[f*g](t) =\n"); pretty (fg, "ascii");
        disp ("\n\nLaTex code (to use in docs):\n");
        disp ("f(t):");disp(latex(f));disp("");
        disp ("f(x-t):");disp(latex(fxt));disp("");
        disp ("g(t):");disp(latex(g));disp("");
        disp ("f*g:");disp(latex(fg));disp("");
    end
    if(opt=="h")
        disp("\n\nTransform properties:") ########TO-DO <-----
        disp("TO-DO")
        disp("\nTransform table:")
        f=dirac(t)
        F=laplace (f,t,s)
        disp("\n")
        f=heaviside(t-n)
        F=laplace (f,t,s)
        f=1
        F=laplace (f,t,s)
        f=t^n
        F=laplace (f,t,s)
        f=exp(-n*t)
        F=laplace (f,t,s)
        f=t*exp(-n*t)
        F=laplace (f,t,s)
        f=cos(m*t)
        F=laplace (f,t,s)
        f=sin(m*t)
        F=laplace (f,t,s)
        f=exp(-n*t)*cos(m*t)
        F=laplace (f,t,s)
        f=exp(-n*t)*sin(m*t)
        F=laplace (f,t,s)
        f=cosh(n*t)
        F=laplace (f,t,s)
        f=sinh(n*t)
        F=laplace (f,t,s)
    end
    if(opt=="e")
        cont = 0;disp ("\n\n\Bye! :D\n");exit;
    end
end