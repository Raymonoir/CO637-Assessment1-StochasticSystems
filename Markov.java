import java.lang.Math;
import java.util.*;

class Markov {

  public static int aproxInfinity = 80000000;

  // Q1 && Q2
  // Get transition probability from s1 to s2 on a grid with numStates number of
  // states
  public static double getTransProb(int s1, int s2, int numStates) {

    int gridSize = (int) Math.sqrt(numStates);
    boolean canMove = possibleMove(s1, s2, gridSize);

    // First check for self trans
    if (s1 == s2) {
      return getSelfTransProb(s1, gridSize);
    } else {

      if (canMove) {
        double pProp = 1.0 / 4.0; // Up down left right
        double pAcc = pAccFunc(s1, s2, numStates);
        double pTrans = pProp * pAcc;
        return pTrans;

      } else {
        return 0.0;
      }

    }
  }

  // Q3
  // Estimate probability of starting in state s1, ending in state 2 in TS number
  // of steps in a discrete markov chain with numStates number of states
  public static double getSejProb(int s1, int s2, int numStates, int TS) {
    return findProbabilities(aproxInfinity, TS, numStates, s2, s1);
  }

  // Q4
  // Get probability of going form state s1 to state s2 using the steady state
  // probabilities provided in ssprob
  public static double getBiasTransProb(int s1, int s2, double[] ssprob) {
    int gridSize = 3;

    boolean canMove = possibleMove(s1, s2, gridSize);

    if (canMove) {

      // If non self transition, use function to do so
      if (s1 != s2) {
        return calcBiasTransProb(s1, s2, ssprob);
      }
      // if self transition, figure out sum of outgoing transitions, all transitions
      // must sum to 1 so self transition prob = 1 - non self transitions probs
      else {
        boolean[] possibleMoves = getPossibleMoves(s1, gridSize);
        double sumTrans = 0;
        // Can go north
        if (possibleMoves[0]) {
          sumTrans += calcBiasTransProb(s1, s1 + gridSize, ssprob);
        }
        // Can go east
        if (possibleMoves[1]) {
          sumTrans += calcBiasTransProb(s1, s1 + 1, ssprob);
        }
        // can go south
        if (possibleMoves[2]) {
          sumTrans += calcBiasTransProb(s1, s1 - gridSize, ssprob);
        }
        // Can go west
        if (possibleMoves[3]) {
          sumTrans += calcBiasTransProb(s1, s1 - 1, ssprob);
        }

        double selfTransProb = 1 - sumTrans;
        return selfTransProb;
      }

    } else {
      return 0.0;
    }

  }

  // Q5
  // In a continous markov chain, using the rates provided what is the probability
  // of moving from state 1 to state 2
  public static double getContTransProb(int s1, int s2, double[] rates) {

    // Self transitions make no sense in continous time markov chain
    if (s1 == s2) {
      return 0;
    } else {
      // Ugly but figure out exactly what node we are moving from and to and calculate
      // probability accordingly using Ki / (Kj + Ki)

      // from state 1
      if (s1 == 1) {
        if (s2 == 2) {
          return rates[0] / (rates[0] + rates[1]);
        } else {
          return rates[1] / (rates[0] + rates[1]);
        }
      }
      // from state 2
      else if (s1 == 2) {
        if (s2 == 1) {
          return rates[2] / (rates[2] + rates[3]);
        } else {
          return rates[3] / (rates[2] + rates[3]);
        }
      }
      // from state 3
      else {
        if (s2 == 1) {
          return rates[4] / (rates[4] + rates[5]);
        } else {
          return rates[5] / (rates[4] + rates[5]);
        }
      }
    }
  }

  // Q6
  // In a continous markov chain, estimate the probability of moving from state 1
  // to state 2 within the time TSC with the provided rates
  public static double getContSejProb(int s1, int s2, double[] rates, double TSC) {
    double t;

    int steps = aproxInfinity;

    int count[] = new int[4];

    for (int i = 0; i < steps; i++) {
      int nextState = s1;
      for (t = 0; t <= TSC; t += getDT(s1, rates)) {
        nextState = contTowerSample(nextState, rates);
      }

      count[nextState]++;
    }
    return (double) count[s2] / (double) steps;
  }

  // ##################################
  // Calculates pAcc using pAcc = Min(1, pi(b) / pi(a)), where pi(a) == ssprob of
  // being in a
  public static double pAccFunc(int s1, int s2, int numStates) {
    double piA = (double) 1 / numStates; // ss prob for grid with numstates = 1/numStates
    double piB = (double) 1 / numStates;
    double piBDivPiA = (double) piB / piA;
    double pAcc = Math.min(1.0, piBDivPiA);

    return pAcc;
  }

  // Returns if the possible move is actually possible ie mvoing from 1 -> 4 is
  // but 1 -> 5 is not in a 3x3 grid
  // 7 8 9
  // 4 5 6
  // 1 2 3
  public static boolean possibleMove(int s1, int s2, int gridSize) {
    boolean canMove = false;
    boolean[] possibleMoves = getPossibleMoves(s1, gridSize);
    // N E S W
    int[] possibleMoveValues = { s1 + gridSize, s1 + 1, s1 - gridSize, s1 - 1 };
    for (int i = 0; i < possibleMoveValues.length; i++) {
      // If the value is in the possible values AND that move is valid (Because going
      // from 3-4 is not a valid move in a 3x3 grid but it will be in the possible
      // values)
      if (possibleMoveValues[i] == s2) {
        if (possibleMoves[i]) {
          canMove = true;
        }
      }
    }
    return canMove;
  }

  // Calculated pAcc similar to that above but uses the provided steady state
  // probabilities
  public static double biasPAccFunc(int s1, int s2, double[] ssProbs) {
    double piA = ssProbs[s1 - 1];
    double piB = ssProbs[s2 - 1];
    double piBDivPiA = (double) piB / piA;
    double pAcc = Math.min(1.0, piBDivPiA);

    return pAcc;
  }

  // Returns a random next state from the current state s1 using the tower
  // sampling method
  public static int towerSample(int s1, int gridSize) {

    boolean[] possibleMoves = getPossibleMoves(s1, gridSize);
    // Probs array in form of {self, N, E, S, W}
    double[] probsArr = { getSelfTransProb(s1, gridSize), 0.25, 0.25, 0.25, 0.25 };

    // if you cannot move in that direction, make the probability of moving in that
    // direction, 0
    for (int i = 1; i < probsArr.length; i++) {
      if (!possibleMoves[i - 1]) {
        probsArr[i] = 0;
      }
    }

    // Random number between 1 and 0
    double r = Math.random();

    double[] cumulProbArray = new double[5];

    cumulProbArray[0] = probsArr[0];
    // Create the cumulative array, using this means if the chance is 0, it is
    // ignored
    for (int i = 1; i < cumulProbArray.length; i++) {
      cumulProbArray[i] = cumulProbArray[i - 1] + probsArr[i];
    }

    // The index of the move I will take in the probs array
    // 0 = self, 1 = N, 2 = E, 3 = S, 4 = W
    int move = 0;

    while (cumulProbArray[move] < r) {
      // If move + 1 is the correct index, set it to that and then break
      if (r <= cumulProbArray[move + 1]) {
        move++; // We add 1 here due to the us picking the upper bound once the index is found
        break;

      } else {
        move++;
      }
    }

    // Self trans
    if (move == 0) {
      return s1;
    }
    // N = current + gridSize
    else if (move == 1) {
      return s1 + gridSize;
    }
    // E = current + 1
    else if (move == 2) {
      return s1 + 1;
    }
    // S = current - gridSize
    else if (move == 3) {
      return s1 - gridSize;
    }
    // W = current - 1
    else {
      return s1 - 1;
    }
  }

  // Returns boolean array in form N E S W if the current node s1 can make that
  // move
  public static boolean[] getPossibleMoves(int s1, int gridSize) {
    // N E S W == up right down left
    boolean[] allMoves = { true, true, true, true };

    if (s1 <= gridSize) {
      // On bottom row, can't move down
      allMoves[2] = false;
    }
    if (s1 % gridSize == 0) {
      // On right wall, cant go right
      allMoves[1] = false;
    }
    if ((s1 - 1) % gridSize == 0) {
      // on left wall, cant go left
      allMoves[3] = false;
    }
    if (s1 > (gridSize - 1) * gridSize) {
      // On top row, cant go up
      allMoves[0] = false;
    }

    return allMoves;
  }

  // This algorithm completes a random walk which is of discrete time 'inner'
  // Once the walk is over, the count of the end state of the walk is increased by
  // 1
  // This is then done 'outer' number of times (hoping to simulate infinity)
  // You can then divide the number of times you landed on a state by 'outer'
  // which will give you a probability of ending on a state in 'inner' time
  // (discrete) - Simulating steady state probabilities
  public static double findProbabilities(int outer, double inner, int k, int probState, int start) {
    int[] count = new int[k + 1];

    for (int i = 0; i < outer; i++) {
      int nextState = start;
      for (int j = 0; j < inner; j++) {
        nextState = towerSample(nextState, (int) Math.sqrt(k));
      }
      count[nextState]++;
    }

    return ((double) count[probState] / (double) outer);
  }

  // Function to calculate the chance of transitioning to self when in a grid
  // 7 8 9
  // 4 5 6
  // 1 2 3
  // Also works for 2x2 as long as you can still only move NESW
  // s1 = Current node
  // gridSize = grid size ie k * k
  public static double getSelfTransProb(int s1, int gridSize) {
    double prob = 0;
    if (s1 <= gridSize) {
      // On bottom row, can't move down
      prob += 0.25;
    }
    if (s1 % gridSize == 0) {
      // On right wall, cant go right
      prob += 0.25;
    }
    if ((s1 - 1) % gridSize == 0) {
      // on left wall, cant go left
      prob += 0.25;
    }
    if (s1 > (gridSize - 1) * gridSize) {
      // On top row, cant go up
      prob += 0.25;
    }

    return prob;
  }

  // Calculates the transition probability from s1 to s2 using the steady state
  // probability provided
  public static double calcBiasTransProb(int s1, int s2, double[] ssprob) {
    double pProp = 1.0 / 4.0; // Up down left right
    double pAcc = biasPAccFunc(s1, s2, ssprob);
    double pTrans = pProp * pAcc;
    return pTrans;
  }

  // Returns a tower sample from state s1 of a continous markov chain with
  // provided rates and only 3 states
  public static int contTowerSample(int s1, double[] rates) {

    int k = 2; // Number of transitions, 2 as they don't have self trans and there are 3 nodes
    Double[] probsArray = new Double[k];

    // Places the correct probabilities in the probs array from sampling
    if (s1 == 1) {
      probsArray[0] = getContTransProb(s1, 2, rates);
      probsArray[1] = getContTransProb(s1, 2, rates);
    } else if (s1 == 2) {
      probsArray[0] = getContTransProb(s1, 1, rates);
      probsArray[1] = getContTransProb(s1, 3, rates);
    } else {
      probsArray[0] = getContTransProb(s1, 1, rates);
      probsArray[1] = getContTransProb(s1, 2, rates);
    }

    double r = Math.random();

    int move = 0;

    // Uses r to select random state depending on transition states
    if (r > probsArray[0]) {
      move = 1;
    }

    // Return correct next state depending on the assigned move
    if (s1 == 1) {
      if (move == 0) {
        return 2;
      } else {
        return 3;
      }
    } else if (s1 == 2) {

      if (move == 0) {
        return 1;
      } else {
        return 3;
      }

    } else {
      if (move == 0) {
        return 1;
      } else {
        return 2;
      }

    }
  }

  // Gets the change in time for a continous markov chain depening on the current
  // state s1 and provided rates, using the function:
  // dt = (-1 / P) * ln(r) where r is a random number between 0 and 1 and P is the
  // sum of the rates out of the state
  public static double getDT(int s1, double[] rates) {
    double transSum = 0;

    if (s1 == 1) {
      transSum += rates[0];
      transSum += rates[1];
    } else if (s1 == 2) {
      transSum += rates[2];
      transSum += rates[3];
    } else {
      transSum += rates[4];
      transSum += rates[5];
    }

    double r = Math.random();
    double dt = (-1.0 / transSum) * Math.log(r);
    return dt;

  }
}
// end class
