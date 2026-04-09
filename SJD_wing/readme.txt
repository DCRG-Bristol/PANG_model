* Black-box models for numerical predictions:
- 'model_G1_1.m' (for test G1.1, i.e., ground static test, clean wing)
- 'model_G1_2.m' (for test G1.2, i.e., ground static test with tip weight)
- 'model_G2.m' (for test G2, i.e., ground vibration test)
- 'model_W1_W2.m' (for tests W1.1, W1.2 and W2, i.e., wind tunnel tests)
- each model is a MATLAB function which takes as inputs the uncertain parameters
* To extract the experimental data for the ground test:
- load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
* For direct comparisons between numerical and experimental data for the ground test:
- in your black-box models, fix the wing root pitch angles to match the 'rootAngl' from 'exprData'
- for visual comparisons, you can run the script 'example_groundTests.m' 
* Additional comments and instructions are given in each function and script.

