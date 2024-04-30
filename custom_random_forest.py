from multiprocessing import Pool

import numpy as np
from numpy import ndarray
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier


class RandomForestClassifierCustom(BaseEstimator):
    """
    Custom Random Forest Classifier.

    Parameters:
        n_estimators (int): The number of trees in the forest.
        max_depth (int): The maximum depth of the tree.
        max_features (int): The number of features to consider when looking for the best split.
        random_state (int): The seed used by the random number generator.

    Attributes:
        classes_ (ndarray): The classes labels.
        n_estimators (int): The number of trees in the forest.
        max_depth (int): The maximum depth of the tree.
        max_features (int): The number of features to consider when looking for the best split.
        random_state (int): The seed used by the random number generator.
        trees (list): List of decision tree classifiers.
        feat_ids_by_tree (list): List of feature indices used by each tree.

    Methods:
        fit(X: ndarray, y: ndarray, n_jobs: int = 1) -> None: Fits the model to the training data.
        predict_proba(X: ndarray, n_jobs: int = 1) -> ndarray: Predicts class probabilities for input data.
        predict(X: ndarray, n_jobs: int) -> ndarray: Predicts class labels for input data.
    """

    def __init__(
        self,
        n_estimators: int = 10,
        max_depth: int = None,
        max_features: int = None,
        random_state: int = 111,
    ) -> None:
        self.classes_ = None
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def fit(self, X: ndarray, y: ndarray, n_jobs: int = 1) -> None:
        """
        Fits the model to the training data.

        Parameters:
            X (ndarray): The training input samples.
            y (ndarray): The target values.
            n_jobs (int): The number of jobs to run in parallel. Defaults to 1.

        Returns:
            None
        """
        self.classes_ = sorted(np.unique(y))

        pool = Pool(processes=n_jobs)
        results = []

        for i in range(self.n_estimators):
            np.random.seed(self.random_state + i)
            features = np.random.choice(X.shape[1], self.max_features, replace=False)
            self.feat_ids_by_tree.append(features)

            bootstrap_idx = np.random.choice(
                X.shape[0],
                size=X.shape[0],
                replace=True,
            )

            X_bagging = X[bootstrap_idx][:, features]
            y_bagging = y[bootstrap_idx]

            if X_bagging.ndim == 1:
                X_bagging = X_bagging.reshape(-1, 1)

            tree_classifier = DecisionTreeClassifier(
                max_depth=self.max_depth,
                max_features=self.max_features,
                random_state=self.random_state,
            )

            results.append(
                pool.apply_async(tree_classifier.fit, (X_bagging, y_bagging))
            )

        pool.close()
        pool.join()

        self.trees = [res.get() for res in results]

    def predict_proba(self, X: ndarray, n_jobs: int = 1) -> ndarray:
        """
        Predicts class probabilities for input data.

        Parameters:
            X (ndarray): The input samples.
            n_jobs (int): The number of jobs to run in parallel. Defaults to 1.

        Returns:
            ndarray: The class probabilities.
        """
        pool = Pool(processes=n_jobs)
        results = []

        for idx, tree in enumerate(self.trees):
            features = self.feat_ids_by_tree[idx]
            X_bagging = X[:, features]

            if X_bagging.ndim == 1:
                X_bagging = X_bagging.reshape(-1, 1)

            results.append(pool.apply_async(tree.predict_proba, (X_bagging,)))

        pool.close()
        pool.join()

        probabilities = [res.get() for res in results]
        return np.mean(probabilities, axis=0)

    def predict(self, X: ndarray, n_jobs: int) -> ndarray:
        """
        Predicts class labels for input data.

        Parameters:
            X (ndarray): The input samples.
            n_jobs (int): The number of jobs to run in parallel.

        Returns:
            ndarray: The predicted class labels.
        """
        probs = self.predict_proba(X, n_jobs=n_jobs)
        predictions = np.argmax(probs, axis=1)

        return predictions
