import sys
import numpy as np
import time
import gtsam
import gtsam.utils.plot as gtsam_plot
import matplotlib.pyplot as plt
from functools import partial
from typing import List, Optional

np.set_printoptions(linewidth=400)

def car():
    # Create noise models
    ODOMETRY_NOISE = gtsam.noiseModel.Diagonal.Sigmas(np.array([0.2, 0.2, 0.1]))
    PRIOR_NOISE = gtsam.noiseModel.Diagonal.Sigmas(np.array([0.3, 0.3, 0.1]))

    graph = gtsam.NonlinearFactorGraph()
    priorMean = gtsam.Pose2(0.0, 0.0, 0.0)
    graph.add(gtsam.PriorFactorPose2(1, priorMean, PRIOR_NOISE))
    odometry = gtsam.Pose2(2.0, 0.0, 0.0)
    variables = 50
    for i in range(1, variables):
        graph.add(gtsam.BetweenFactorPose2(i, i + 1, odometry, ODOMETRY_NOISE))

    print("\nFactor Graph:\n{}".format(graph))
    initial = gtsam.Values()
    for i in range(1, variables + 1):
        initial.insert(i, gtsam.Pose2((i - 1) * 2 + np.random.normal(0, 0.2, 1), np.random.normal(0, 0.2, 1), np.random.normal(0, 0.1, 1)))
    print("\nInitial Estimate:\n{}".format(initial))
    # params = gtsam.LevenbergMarquardtParams()
    params = gtsam.GaussNewtonParams()
    params.setRelativeErrorTol(1e-5)
    params.setMaxIterations(100)
    params.setVerbosity("ERROR".encode("utf-8"))
    # optimizer = gtsam.LevenbergMarquardtOptimizer(graph, initial, params)

        
    optimizer = gtsam.GaussNewtonOptimizer(graph, initial, params)
    result = optimizer.optimize()



    print("\nFinal Result:\n{}".format(result))
    marginals = gtsam.Marginals(graph, result)
    for i in range(1, variables + 1):
        print("X{} covariance:\n{}\n".format(i,
                                             marginals.marginalCovariance(i)))
        
    for i in range(1, variables + 1):
        gtsam_plot.plot_pose2(0, initial.atPose2(i), 0.5)
    plt.axis('equal')

    for i in range(1, variables + 1):
        # gtsam_plot.plot_pose2(1, result.atPose2(i), 0.5, marginals.marginalCovariance(i))
        gtsam_plot.plot_pose2(1, result.atPose2(i), 0.5)
    plt.axis('equal')
    plt.show()

def error_point(measurement: np.ndarray, this: gtsam.CustomFactor,
               values: gtsam.Values,
               jacobians: Optional[List[np.ndarray]]) -> float:
    """Point (only position, no orientation) Factor error function
    :param measurement: Point measurement, to be filled with `partial`
    :param this: gtsam.CustomFactor handle
    :param values: gtsam.Values
    :param jacobians: Optional list of Jacobians
    :return: the unwhitened error
    """
    key1 = this.keys()[0]
    key2 = this.keys()[1]
    pos1, pos2 = values.atVector(key1), values.atVector(key2)
    error =   measurement - (pos2 - pos1)
    if jacobians is not None:
        jacobians[0] = np.eye(2)
        jacobians[1] = -np.eye(2)

    return error
    
def point():
    ODOMETRY_NOISE = gtsam.noiseModel.Diagonal.Sigmas(np.array([0.2, 0.2]))
    PRIOR_NOISE = gtsam.noiseModel.Diagonal.Sigmas(np.array([0.3, 0.3]))
    latency = []
    graph = gtsam.NonlinearFactorGraph()

    priorMean = gtsam.Point2(0.0, 0.0)
    graph.add(gtsam.PriorFactorPoint2(1, priorMean, PRIOR_NOISE))

    # odometry = gtsam.Point2(2.0, 0.0)
    variables = 50

    initial = gtsam.Values()
    for i in range(1, variables + 1):
        initial.insert(i, gtsam.Point2((i - 1) * 2 + np.random.normal(0, 0.2), np.random.normal(0, 0.2)))
    print("\nInitial Estimate:\n{}".format(initial))

    for i in range(1, variables):
        odof = gtsam.CustomFactor(ODOMETRY_NOISE, [i, i + 1], partial(error_point, np.array([2.0, 0.0])))
        graph.add(odof)
    print("\nFactor Graph:\n{}".format(graph))

    # params = gtsam.LevenbergMarquardtParams()
    params = gtsam.GaussNewtonParams()
    params.setRelativeErrorTol(1e-5)
    params.setMaxIterations(100)
    params.setVerbosity("ERROR".encode("utf-8"))
    # optimizer = gtsam.LevenbergMarquardtOptimizer(graph, initial, params)
 

    optimizer = gtsam.GaussNewtonOptimizer(graph, initial, params)
    result = optimizer.optimize()



    print("\nFinal Result:\n{}".format(result))
    # marginals = gtsam.Marginals(graph, result)
    # for i in range(1, variables + 1):
    #     print("X{} covariance:\n{}\n".format(i, marginals.marginalCovariance(i)))

    initial_x = []
    initial_y = []
    result_x = []
    result_y = []

    for i in range(1, variables + 1):
        initial_x.append(initial.atVector(i)[0])
        initial_y.append(initial.atVector(i)[1])
        result_x.append(result.atVector(i)[0])
        result_y.append(result.atVector(i)[1])

    plt.subplot(2,1,1)
    plt.ylim(-1, 1)
    plt.plot(initial_x, initial_y, marker='o')
    plt.subplot(2,1,2)
    plt.ylim(-1, 1)
    plt.plot(result_x, result_y, marker='*')

    plt.show()


    
        
    # for i in range(1, variables + 1):
    #     p = gtsam.Point3()
    #     p[0:2] = initial.atPoint2(i)
    #     p[2] = 0
    #     print(p)
    #     gtsam_plot.plot_point3(0, p, 'bo')
    # # plt.axis('equal')

    # for i in range(1, variables + 1):
    #     p = gtsam.Point3()
    #     p[0:2] = result.atPoint2(i)
    #     p[2] = 0
    #     print(p)
    #     gtsam_plot.plot_point3(1, p, 'bo')
    # # plt.axis('equal')
    # plt.show()


if __name__ == "__main__":
    
    car()
    point()