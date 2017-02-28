function [integrated_out,NPTS_out,RELERR_out,IFAIL_out]=D01AHF(A,B,EPR,F,NL,IFAIL)
IFAIL_out=0; %Is not used, but later will be nice to use
RELERR_out=0; %Is not used, but later will be nice to use
% [integrated_out,fcnt] = quad(F,A,B,EPR);
% [integrated_out,fcnt] = quad(F,A,B); %this gives the best result but is slow
% [integrated_out,fcnt] = quadgk(F,A,B,'AbsTol',1e-10,'RelTol', 0, 'MaxIntervalCount',20000 ); %the faster
[integrated_out,fcnt] = quadgk(F,A,B); %the faster
% [integrated_out,fcnt] = quadv(F,A,B); %slow also
% [integrated_out]=rombint(F,A,B,15);
% [integrated_out]=quadgr(F,A,B,1e-09); %nada mal
% fcnt=1;
NPTS_out=fcnt;
    
%THIS FUNCTION ROUTINE PERFORMS AUTOMATIC INTEGRATION OVER A FINITE INTERVAL 
%USING THE BASIC INTEGRATION ALGORITHMS 
%     INPUT ARGUMENTS
%     ----- ----------
%     A,B     -  LOWER AND UPPER INTEGRATION LIMITS.
%     EPR     -  REQUIRED RELATIVE ACCURACY.
%     NL      -  APPROXIMATE LIMIT ON NUMBER OF INTEGRAND
%                EVALUATIONS. IF SET NEGATIVE OR ZERO THE
%                DEFAULT IS 10000.
%     F       -  THE USER NAMED AND PREPARED FUNCTION  F(X)
%                GIVES THE VALUE OF THE INTEGRAND AT X.
%     IFAIL      INTEGER VARIABLE
%             - 0  FOR HARD FAIL REPORT
%             - 1  FOR SOFT FAIL REPORT
%    OUTPUT ARGUMENTS
%     ------ ----------
%     NPTS    -  NUMBER OF INTEGRAND EVALUATIONS USED IN OBTAINING
%                THE RESULT.
%     RELERR  -  ROUGH ESTIMATE OF RELATIVE ACCURACY ACHIEVED.
%     IFAIL   -  VALUE INDICATES THE OUTCOME OF THE INTEGRATION -
%                IFAIL  = 0  CONVERGED
%                IFAIL  = 1  INTEGRAND EVALUATIONS EXCEEDED  NL.
%                            THE RESULT WAS OBTAINED BY CONTINUING
%                            BUT IGNORING ANY NEED TO SUBDIVIDE.
%                            RESULT LIKELY TO BE INACCURATE.
%                IFAIL  = 2  DURING THE SUBDIVISION PROCESS
%                            THE STACK BECAME FULL
%                            (PRESENTLY SET TO HOLD 20
%                            LEVELS OF INFORMATION.  MAY BE
%                            INCREASED BY  ALTERING  ISMAX
%                            AND THE DIMENSIONS OF STACK
%                            AND ISTACK). RESULT IS
%                            OBTAINED BY CONTINUING BUT
%                            IGNORING CONVERGENCE FAILURES
%                            ON INTERVALS  WHICH CANNOT BE
%                            ACCOMMODATED ON THE STACKS.
%                            RESULT LIKELY TO BE
%                            INACCURATE.
%                IFAIL  = 3  INVALID ACCURACY REQUEST.