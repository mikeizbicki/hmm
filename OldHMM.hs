{-# LANGUAGE ParallelListComp #-}

module OldHMM
    (Prob, HMM(..), train, bestSequence, sequenceProb)
    where

import qualified Data.Map as M
import Data.List (sort, groupBy, maximumBy, foldl')
import Data.Maybe (fromMaybe, fromJust)
import Data.Ord (comparing)
import Data.Function (on)
import Control.Monad
import Data.Number.LogFloat

type Prob = LogFloat

-- | The type of Hidden Markov Models.
data HMM state observation = HMM [state] [Prob] [[Prob]] (observation -> [Prob])

instance (Show state, Show observation) => Show (HMM state observation) where
    show (HMM states probs tpm _) = "HMM " ++ show states ++ " "
                                           ++ show probs ++ " " ++ show tpm ++ " <func>"

-- | Perform a single step in the Viterbi algorithm.
--  
--   Takes a list of path probabilities, and an observation, and returns the updated
--   list of (surviving) paths with probabilities.
viterbi ::     HMM state observation
            -> [(Prob, [state])]
            -> observation
            -> [(Prob, [state])]
viterbi (HMM states _ state_transitions observations) prev x =
    deepSeq prev `seq`
    [maximumBy (comparing fst)
            [(transition_prob * prev_prob * observation_prob,
               new_state:path)
                    | transition_prob <- transition_probs
                    | (prev_prob, path) <- prev
                    | observation_prob <- observation_probs]
        | transition_probs <- state_transitions
        | new_state <- states]
    where
        observation_probs = observations x
        deepSeq ((x, y:ys):xs) = x `seq` y `seq` (deepSeq xs)
        deepSeq ((x, _):xs) = x `seq` (deepSeq xs)
        deepSeq [] = []

-- | The initial value for the Viterbi algorithm
viterbi_init :: HMM state observation -> [(Prob, [state])]
viterbi_init (HMM states state_probs _ _) = zip state_probs (map (:[]) states)

-- | Perform a single step of the forward algorithm
-- 
--   Each item in the input and output list is the probability that the system
--   ended in the respective state.
forward ::     HMM state observation
            -> [Prob]
            -> observation
            -> [Prob]
forward (HMM _ _ state_transitions observations) prev x =
    last prev `seq`
    [sum [transition_prob * prev_prob * observation_prob
                | transition_prob <- transition_probs
                | prev_prob <- prev
                | observation_prob <- observation_probs]
        | transition_probs <- state_transitions]
    where
        observation_probs = observations x

-- | The initial value for the forward algorithm
forward_init :: HMM state observation -> [Prob]
forward_init (HMM _ state_probs _ _) = state_probs

learn_states :: (Ord state, Fractional prob) => [(observation, state)] -> M.Map state prob
learn_states xs = histogram $ map snd xs

learn_transitions :: (Ord state, Fractional prob) => [(observation, state)] -> M.Map (state, state) prob
learn_transitions xs = let xs' = map snd xs in
                        histogram $ zip xs' (tail xs')

learn_observations ::  (Ord state, Ord observation, Fractional prob) =>
                       M.Map state prob
                    -> [(observation, state)]
                    -> M.Map (observation, state) prob
learn_observations state_prob = M.mapWithKey (\ (observation, state) prob -> prob / (fromJust $ M.lookup state state_prob))
                            . histogram

histogram :: (Ord a, Fractional prob) => [a] -> M.Map a prob
histogram xs = let hist = foldl' (flip $ flip (M.insertWith (+)) 1) M.empty xs in
                M.map (/ M.fold (+) 0 hist) hist

-- | Calculate the parameters of an HMM from a list of observations
--   and the corresponding states.
train :: (Ord observation, Ord state) =>
            [(observation, state)]
         -> HMM state observation
train sample = model
    where
        states = learn_states sample
        state_list = M.keys states
        
        transitions = learn_transitions sample
        trans_prob_mtx = [[fromMaybe 1e-10 $ M.lookup (old_state, new_state) transitions
                                | old_state <- state_list]
                                | new_state <- state_list]

        observations = learn_observations states sample
        observation_probs = fromMaybe (fill state_list []) . (flip M.lookup $
                            M.fromList $ map (\ (e, xs) -> (e, fill state_list xs)) $
                                map (\ xs -> (fst $ head xs, map snd xs)) $
                                groupBy     ((==) `on` fst)
                                            [(observation, (state, prob))
                                                | ((observation, state), prob) <- M.toAscList observations])

        initial = map (\ state -> (fromJust $ M.lookup state states, [state])) state_list

        model = HMM state_list (fill state_list $ M.toAscList states) trans_prob_mtx observation_probs

        fill :: Eq state => [state] -> [(state, Prob)] -> [Prob]
        fill states [] = map (const 1e-10) states
        fill (s:states) xs@((s', p):xs') = if s /= s' then
                                            1e-10 : fill states xs
                                           else
                                            p : fill states xs'

-- | Calculate the most likely sequence of states for a given sequence of observations
--   using Viterbi's algorithm
bestSequence :: (Ord observation) => HMM state observation -> [observation] -> [state]
bestSequence hmm = (reverse . tail . snd . (maximumBy (comparing fst))) . (foldl' (viterbi hmm) (viterbi_init hmm))

-- | Calculate the probability of a given sequence of observations
--    using the forward algorithm.
sequenceProb :: (Ord observation) => HMM state observation -> [observation] -> Prob
sequenceProb hmm = sum . (foldl' (forward hmm) (forward_init hmm))