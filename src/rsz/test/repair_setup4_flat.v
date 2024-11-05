module reg1 (clk);
 input clk;

 wire r1q;
 wire u1z;
 wire u2z;
 wire u3z;
 wire u4z;
 wire u5z;

 DFF_X1 r1 (.CK(clk),
    .Q(r1q));
 BUF_X1 u1 (.A(r1q),
    .Z(u1z));
 BUF_X1 u2 (.A(u1z),
    .Z(u2z));
 BUF_X1 u3 (.A(u2z),
    .Z(u3z));
 BUF_X1 u4 (.A(u3z),
    .Z(u4z));
 BUF_X1 u5 (.A(u4z),
    .Z(u5z));
 DFF_X1 r2 (.D(u5z),
    .CK(clk));
 DFF_X1 r3 (.D(r1q),
    .CK(clk));
 DFF_X1 r4 (.D(r1q),
    .CK(clk));
 DFF_X1 r5 (.D(r1q),
    .CK(clk));
 DFF_X1 r6 (.D(r1q),
    .CK(clk));
 DFF_X1 r7 (.D(r1q),
    .CK(clk));
 DFF_X1 r8 (.D(r1q));
 DFF_X1 r9 (.D(r1q));
 DFF_X1 r10 (.D(r1q));
 DFF_X1 r11 (.D(r1q));
 DFF_X1 r12 (.D(r1q));
endmodule
