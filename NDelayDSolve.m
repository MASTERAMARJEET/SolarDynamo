(* :Title:NDelayDSolve *)
(* :Author: Allan Hayes, hay@haystack.demon.co.uk *)
(* :Summary:
    NDelayDSolve is a package containing a single 
    function, NDelayDSolve, for the numerical solution 
    of delay-differential equations.
*)
(* :Context: haypacks`NumericalMath`NDelayDSolve` *)
(* :Package Version: 1.0 *)
(* :Copyright: Copyright 1997 Allan Hayes.
Permission to duplicate this work is granted provided 
that the copyright notice is included. 
The author shall retain the copyrght of this work*)

(* :History:        
   Version 0.1 by Allan Hayes, December 1997 
   Version 1.0 including changes by Mats Jirstrand to deal with new features of Mathematica 5.0, January 2004  
*)
(* :Warning: 
*)
(* :Keywords: Delay, DifferentialEquation, NDSolve *)
(* :Mathematica Version: 3.01 *)
(* :Limitation: *)

(**Begin Package**)

BeginPackage[
    "haypacks`NumericalMath`NDelayDSolve`",
    "Utilities`FilterOptions`"];

Unprotect["`*"];ClearAll["`*"];

(**usage messages**)

NDelayDSolve::usage = "\n
NDelayDSolve[eqns, initsoln, {t,tmin,tmax}]\n 
gives a numerical solution, as a list of 
interpolating functions, for a list of ordinary and
delayed differential equations. The input pattern, 
which is slightly different from that for NDSolve, 
is shown by the following examples:\n
NDelayDSolve[
    {x''[t] == 5*(t-.9)^2 -3(x')[t-1.9]},
    {x->(#1&)},{t,1,10}]
\n
NDelayDSolve[
    { x''[t] == (1+y[t-1.1])-5x'[t-1.9],
      y''[t] == x[t-2.7]+ 6y'[t-3.5]
    },
    {x->(1&),y->(0&)},{t,0,8}
]
\n
The initial functions, given by rules like x ->  function (eg 1& or Sin), \
must allow for the maximum delay (in the first example, over -0.9 to 1 at \
least)

Equations without delay have to be input in the same 
way: 
\n
NDelayDSolve[{x'[t] ==12Sin[12 t] + x[t]}, 
   {x->Cos}, {t, 0,1}]
\n gives the same solution as
\n
NDSolve[{x'[t] ==12Sin[12 t] + x[t],x[0]==Cos[0]}, 
   x, {t,0,1}]
\n but with one list level only.\n
Options are passed on to NDSolve and FunctionInterpolation.";

Begin["`Private`"];

Options[NDelayDSolve]=
    Union[Options[NDSolve],Options[FunctionInterpolation]];

NDelayDSolve[eqns_,initsoln_, {t_,tmin_,tmax_},
opts___?OptionQ]:=
Module[
    {vars,altvars,soln,delmin,initexprs,func,
    neweqns,initvals,b,tsup},

    (** set up **)
        
    vars= First/@initsoln;
    altvars = Alternatives@@vars;
    soln[s_/;s<= tmin] = initsoln;
    delmin =Min[Cases[eqns,v_[t+d_]:>-d ,Infinity]];
    
    (* find the expressions that will have to be given 
    initial values*)
        
    initexprs = 
        Cases[eqns,
            Derivative[n_][v:altvars][_]:>
                Table[Derivative[k][v],{k,0,n-1}],
            Infinity
        ]//Flatten//Union;

    (* replace expressions like Derivative[n][x][s] 
    by the following forms that only evaluate when
    s is numeric*)
    
    (* ToExpression inserted in func definition,
    Mats Jirstrand 2004-01-29*)    
    func[u_,n_,s_?NumericQ]:=Derivative[n][ToExpression[u]/.soln[s]][s];

    (** make the replacemnts **)
    
    (* ToString inserted in calls to func,
    Mats Jirstrand 2004-01-29*)      
    neweqns:=
        (eqns/.
            {
                (Derivative[n_][u:altvars][t+d_]):> 
                    func[ToString[u],n,t+d],
                (u:altvars)[t+d_]:> 
                    func[ToString[u],0,t+d]
            }
        );
    
    (** calculate the initial values **)
    
    initvals:=
        Thread[
            Through/@
                (
                initexprs[b]==(initexprs/.soln[b])[b]
                )
        ];
                
    (** solve the DE **)
        
    b=tmin;
    While[b==tmin||b<tmax&&tsup==b,
        soln[s_/;Evaluate[b<s<=b+delmin]]= 
            NDSolve[
                Join[neweqns,initvals],
                vars,
                {t,b,Min[b+delmin,tmax]},
                FilterOptions[NDSolve,opts]
            ]//First; 
                        
        b+=delmin;
        tsup = (First[vars]/.soln[b])[[1,1,2]]
    ];  
    
    (** form single interpolating function **)

    (* ToString inserted in calls to func,
    Mats Jirstrand 2004-01-29*)           
    Cases[vars, 
        v_ :> 
            (v-> 
                FunctionInterpolation[
                    func[ToString[v],0,t],{t,tmin,tsup},
                    FilterOptions[
                        FunctionInterpolation,opts
                    ]//Evaluate
                                  
                ]
            )
    ]
];

End[];
Protect["`*"];
EndPackage[];
