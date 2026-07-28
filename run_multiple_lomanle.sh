echo "------------------------------------------------------------------------------------------------"
echo "Running Lomanle 500 points, k=30, iter=10"
./run_lomanle_tests.sh results/data/2d/circular_arc_noise_high.csv 30 1 3.0 0.3 0.10 0.90 2.5 10 0.01 2

echo "------------------------------------------------------------------------------------------------"
echo "Running Lomanle 500 points, k=30, iter=50"
./run_lomanle_tests.sh results/data/2d/circular_arc_noise_high.csv 30 1 3.0 0.3 0.10 0.90 2.5 50 0.01 2

echo "------------------------------------------------------------------------------------------------"
echo "Running Lomanle 500 points, k=80, iter=10"
./run_lomanle_tests.sh results/data/2d/circular_arc_noise_high.csv 80 1 3.0 0.3 0.10 0.90 2.5 10 0.01 2

echo "------------------------------------------------------------------------------------------------"
echo "Running Lomanle 500 points, k=80, iter=50"
./run_lomanle_tests.sh results/data/2d/circular_arc_noise_high.csv 80 1 3.0 0.3 0.10 0.90 2.5 50 0.01 2

echo "------------------------------------------------------------------------------------------------"
echo "Running Lomanle 5000 points, k=300, iter=10"
./run_lomanle_tests.sh results/data/2d/5000/circular_arc_2d_noise_high.csv 300 1 3.0 0.3 0.10 0.90 2.5 10 0.01 2

echo "------------------------------------------------------------------------------------------------"
echo "Running Lomanle 5000 points, k=300, iter=50"
./run_lomanle_tests.sh results/data/2d/5000/circular_arc_2d_noise_high.csv 300 1 3.0 0.3 0.10 0.90 2.5 50 0.01 2

echo "------------------------------------------------------------------------------------------------"
echo "Running Lomanle 5000 points, k=800, iter=10"
./run_lomanle_tests.sh results/data/2d/5000/circular_arc_2d_noise_high.csv 800 1 3.0 0.3 0.10 0.90 2.5 10 0.01 2

echo "------------------------------------------------------------------------------------------------"
echo "Running Lomanle 5000 points, k=800, iter=50"
./run_lomanle_tests.sh results/data/2d/5000/circular_arc_2d_noise_high.csv 800 1 3.0 0.3 0.10 0.90 2.5 50 0.01 2

echo "------------------------------------------------------------------------------------------------"
echo "Running Lomanle 50000 points, k=3000, iter=10"
./run_lomanle_tests.sh results/data/2d/50000/circular_arc_2d_noise_high.csv 3000 1 3.0 0.3 0.10 0.90 2.5 10 0.01 2

echo "------------------------------------------------------------------------------------------------"
echo "Running Lomanle 50000 points, k=3000, iter=50"
./run_lomanle_tests.sh results/data/2d/50000/circular_arc_2d_noise_high.csv 3000 1 3.0 0.3 0.10 0.90 2.5 50 0.01 2

echo "------------------------------------------------------------------------------------------------"
echo "Running Lomanle 50000 points, k=8000, iter=10"
./run_lomanle_tests.sh results/data/2d/50000/circular_arc_2d_noise_high.csv 8000 1 3.0 0.3 0.10 0.90 2.5 10 0.01 2

echo "------------------------------------------------------------------------------------------------"
echo "Running Lomanle 50000 points, k=8000, iter=50"
./run_lomanle_tests.sh results/data/2d/50000/circular_arc_2d_noise_high.csv 8000 1 3.0 0.3 0.10 0.90 2.5 50 0.01 2