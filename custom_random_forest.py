import os
import random
from multiprocessing import Pool

import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier


class RandomForestClassifierCustom(BaseEstimator):
    def __init__(
        self,
        n_estimators=10,
        max_depth=None,
        max_features=None,
        random_state=os.getenv("SEED"),
    ):
        self.classes_ = None
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def fit(self, X, y, n_jobs=1):
        self.classes_ = sorted(np.unique(y))

        pool = Pool(processes=n_jobs)
        results = []

        for i in range(self.n_estimators):
            random.seed(self.random_state + i)
            features = np.random.choice(X.shape[1], self.max_features, replace=False)
            self.feat_ids_by_tree.append(features)

            bootstrap_idx = np.random.choice(X.shape[0], size=X.shape[0], replace=True)

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

    def predict_proba(self, X, n_jobs=1):

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

    def predict(self, X, n_jobs):
        probs = self.predict_proba(X, n_jobs=n_jobs)
        predictions = np.argmax(probs, axis=1)

        return predictions
