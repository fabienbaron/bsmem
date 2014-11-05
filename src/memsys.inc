      REAL RCOM
      INTEGER ICOM,MRCOM,MICOM
      PARAMETER (MRCOM=40,MICOM=100)
      COMMON /MEMCOM/ RCOM(MRCOM),ICOM(MICOM)

* NSIZE  is size of Omega(alpha) table
* LMAX   is maximum number of conjugate gradient iterations
* NINTS  is number of INTEGERs to be preserved in save/restore - useless now
* NREALS is number of REALs to be preserved in save/restore - useless now
* Default values when taking over bsmem : 8, 5, 30, 33
      INTEGER NSIZE,LMAX,NINTS,NREALS
      PARAMETER (NSIZE=8,LMAX=30,NINTS=30,NREALS=33)

      INTEGER KMETH0,KMETH1,KMETH2,KMETH3,KNRAND,KISEED,KNTRNS,
     *        KISTAT,KNTABL,KIOFF,KMORE,KLNGTH,KKC,KKF
      PARAMETER (KMETH0= 1,KMETH1= 2,KMETH2= 3,KMETH3= 4,KNRAND= 5,
     *           KISEED= 6,KNTRNS= 7,KISTAT= 8,KNTABL= 9,KIOFF =10,
     *           KKC   =11,KKF   =51,KMORE =91,KLNGTH=92)
      INTEGER KAIM,KRATE,KDEF,KACC,KUTOL,KTOL,KALPHA,KSCALE,KHHIGH,
     *        KXTABL,KYTABL,KVTABL
      PARAMETER (KAIM  = 1,KRATE = 2,KDEF  = 3,KACC  = 4,KUTOL = 5,
     *           KTOL  = 6,KALPHA= 7,KXTABL= 8,KYTABL=16,KVTABL=24,
     *           KSCALE=32,KHHIGH=33)

