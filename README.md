# Read Me

@last update: 4.12

## HomeWork 4.13-4.20

1. EPNP
According to the paper:[EPNP](https://icwww.epfl.ch/~lepetit/papers/lepetit_ijcv08.pdf),complete the EPNP functions. 
> Lepetit V, Moreno-Noguer F, Fua P. Epnp: An accurate o (n) solution to the pnp problem[J]. International journal of computer vision, 2009, 81(2): 155.

2. Including

- Compiling with Cmake tool chain.
- Available Libray: Eigen , openCV
- Unit Testing: google-test

3. to be Submitted

- Project EPNP
- summary report(PPT) 

✅ 环境配置 cmake,opencv343,eigen
✅ feature extraction, matching
✅ opencv solver,eigen solver
todo: gtest, PPT summary
![img](https://wx3.sinaimg.cn/mw690/c7716318ly1g2172ibxq1j20ze0dx7wh.jpg)


> result: 从打印结果R，T看，两种解法结果基本相同
-- Max dist : 95.000000 
-- Min dist : 7.000000 
一共找到了81组匹配点
3d-2d pairs: 77

---using opencv solvePNP---

R1: [0.9978654469053291, -0.05171629787598182, 0.039874483149384;
 0.0505533717792531, 0.9982813811698801, 0.02964187260119196;
 -0.04133892202484731, -0.02756281087914287, 0.9987649297919227]
t1:[-0.1266609334454444;
 -0.0111141371771197;
 0.05673412814657741]

 ---using egien to solve EPNP---

R1:   0.997854 -0.0518631  0.0399726
 0.0506961   0.998273  0.0296762
-0.0414427 -0.0275861    0.99876
t1: -0.127002
-0.0109942
 0.0567606

