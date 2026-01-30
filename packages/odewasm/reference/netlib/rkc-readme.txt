I'd like to point out an inconsistency in the file https://netlib.org/ode/rkc.f.

Details in the email thread below. 


The fix needed is on line 38, the documentation stating

    If INFO(2) = 0,

should be replaced with

    If INFO(2) = 1,

The memory requirements are stated correctly later in the listing for the argument `WORK(*)`.

I noticed a documentation error in the RKC solver available from Netlib (https://netlib.org/ode/rkc.f).


The documentation states that

c  FIRST CALL TO RKC  
c
c  You must provide storage in your calling program for the arrays in the 
c  call list -- Y(NEQN), INFO(4), WORK(8+5*NEQN).  If INFO(2) = 0, you can
c  reduce the storage for the work array to WORK(8+4*NEQN).

However later in the individual explanation of the arguments, it says

c  WORK(*):  Work array.  Double precision array of length at least 
c            8 + 5*NEQN if INFO(2) = 0 and otherwise, 8 + 4*NEQN.


I have studied the program closely and also using runtime-checking with the gfortran compiler.

The correct dimensions are stated in the section for `WORK(*)`.


I would advise making the small change in the "FIRST CALL TO RKC" section

c  FIRST CALL TO RKC  
c
c  You must provide storage in your calling program for the arrays in the 
c  call list -- Y(NEQN), INFO(4), WORK(8+5*NEQN).  If INFO(2) = 1, ...

Best regards,
Ivan Pribec
