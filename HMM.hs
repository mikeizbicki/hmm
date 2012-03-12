module HMM 
    ( HMM(..), rnf
    , forward
    , backward
    , baumWelch
    , alpha, beta
    )
    where

import Debug.Trace
import Data.Array
import Data.List
import Data.Number.LogFloat
import qualified Data.MemoCombinators as Memo
import Control.DeepSeq
import Control.Parallel

type Prob = LogFloat

   -- | The data type for our HMM

data -- (Eq eventType, Eq stateType, Show eventType, Show stateType) =>
     HMM stateType eventType = HMM { states :: [stateType]
                                   , events :: [eventType]
                                   , initProbs :: (stateType -> Prob)
                                   , transMatrix :: (stateType -> stateType -> Prob)
                                   , outMatrix :: (stateType -> eventType -> Prob)
                                   }

instance NFData (HMM stateType eventType) where
     rnf a = a `seq` ()

instance (Show state, Show observation) => Show (HMM state observation) where
    show hmm = "HMM" ++ " states=" ++ (show $ states hmm) 
                     ++ " events=" ++ (show $ events hmm) 
                     ++ " initProbs=" ++ (show [(s,initProbs hmm s) | s <- states hmm])
                     ++ " transMatrix=" ++ (show [(s1,s2,transMatrix hmm s1 s2) | s1 <- states hmm, s2 <- states hmm])
                     ++ " outMatrix=" ++ (show [(s,e,outMatrix hmm s e) | s <- states hmm, e <- events hmm])

stateIndex :: (Show stateType, Show eventType, Eq stateType) => HMM stateType eventType -> stateType -> Int
stateIndex hmm state = case elemIndex state $ states hmm of 
                            Nothing -> seq (error "stateIndex: Index "++show state++" not in HMM "++show hmm) 0
                            Just x -> x

eventIndex :: (Show stateType, Show eventType, Eq eventType) => HMM stateType eventType -> eventType -> Int
eventIndex hmm event = case elemIndex event $ events hmm of 
                            Nothing -> seq (error "stateIndex: Index "++show event++" not in HMM "++show hmm) 0
                            Just x -> x

   -- | forward algorithm
      
forward :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> [eventType] -> Prob
forward hmm obs = forwardArray hmm (listArray (1,bT) obs)
    where
          bT = length obs
                               
forwardArray :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> Array Int eventType -> Prob
forwardArray hmm obs = sum [alpha hmm obs bT state | state <- states hmm]
    where
          bT = snd $ bounds obs
                                                         
alpha :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType 
                                                                      -> Array Int eventType 
                                                                      -> Int 
                                                                      -> stateType 
                                                                      -> Prob
alpha hmm obs = memo_alpha
    where memo_alpha t state = memo_alpha2 t (stateIndex hmm state)
          memo_alpha2 = (Memo.memo2 Memo.integral Memo.integral memo_alpha3)
          memo_alpha3 t' state'
            | t' == 1       = -- trace ("memo_alpha' t'="++show t'++", state'="++show state') $ 
                              (outMatrix hmm (states hmm !! state') $ obs!t')*(initProbs hmm $ states hmm !! state')
            | otherwise     = -- trace ("memo_alpha' t'="++show t'++", state'="++show state') $ 
                              (outMatrix hmm (states hmm !! state') $ obs!t')*(sum [(memo_alpha (t'-1) state2)*(transMatrix hmm state2 (states hmm !! state')) | state2 <- states hmm])

   
   -- | backwards algorithm
   
backward :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> [eventType] -> Prob
backward hmm obs = backwardArray hmm $ listArray (1,length obs) obs
    
backwardArray :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> Array Int eventType -> Prob
backwardArray hmm obs = backwardArray' hmm obs
    where 
          backwardArray' hmm obs = sum [(initProbs hmm state)
                                       *(outMatrix hmm state $ obs!1)
                                       *(beta hmm obs 1 state)
                                       | state <- states hmm
                                       ]
    
beta :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType 
                                                                      -> Array Int eventType 
                                                                      -> Int 
                                                                      -> stateType 
                                                                      -> Prob
beta hmm obs = memo_beta
    where bT = snd $ bounds obs
          memo_beta t state = memo_beta2 t (stateIndex hmm state)
          memo_beta2 = (Memo.memo2 Memo.integral Memo.integral memo_beta3)
          memo_beta3 t' state'
            | t' == bT       = -- trace ("memo_alpha' t'="++show t'++", state'="++show state') $ 
                              1
            | otherwise     = -- trace ("memo_alpha' t'="++show t'++", state'="++show state') $ 
                              sum [(transMatrix hmm (states hmm !! state') state2)
                                  *(outMatrix hmm state2 $ obs!(t'+1))
                                  *(memo_beta (t'+1) state2) 
                                  | state2 <- states hmm
                                  ]


   -- | Baum-Welch
   
{-gammaArray :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType
                                                                           -> Array Int eventType
                                                                           -> Int
                                                                           -> stateType
                                                                           -> Prob-}
   
   -- xi i j = P(state (t-1) == i && state (t) == j | obs, lambda)
   
-- xiArray :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType 
--                                                                         -> Array Int eventType 
--                                                                         -> Int 
--                                                                         -> stateType 
--                                                                         -> stateType 
--                                                                         -> Prob


baumWelch :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> Array Int eventType -> Int -> HMM stateType eventType
baumWelch hmm obs count
    | count == 0    = hmm
    | otherwise     = itr `seq` baumWelch itr obs (count-1)
        where itr = baumWelchItr hmm obs

baumWelchItr :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> Array Int eventType -> HMM stateType eventType
baumWelchItr hmm obs = --par newInitProbs $ par newTransMatrix $ par newOutMatrix 
                       HMM { states = states hmm
                           , events = events hmm
                           , initProbs = newInitProbs
                           , transMatrix = {-transMatrix hmm-} newTransMatrix
                           , outMatrix = {-outMatrix hmm-} newOutMatrix
                           }
    where bT = snd $ bounds obs
          newInitProbs state = gamma 1 state
          newTransMatrix state1 state2 = sum [xi t state1 state2 | t <- [2..(snd $ bounds obs)]]
                                        /sum [gamma t state1 | t <- [2..(snd $ bounds obs)]]
          newOutMatrix state event = sum [if (obs!t == event) 
                                             then gamma t state 
                                             else 0
                                         | t <- [2..(snd $ bounds obs)]
                                         ]
                                    /sum [gamma t state | t <- [2..(snd $ bounds obs)]]
                                    
          -- Greek functions, included here for memoization
          xi t state1 state2 = (memo_alpha (t-1) state1)
                             *(transMatrix hmm state1 state2)
                             *(outMatrix hmm state2 $ obs!t)
                             *(memo_beta t state2)
                             /backwardArrayVar -- (backwardArray hmm obs)
          
          gamma t state = (memo_alpha t state)
                        *(memo_beta t state)
                        /backwardArrayVar

          backwardArrayVar = (backwardArray hmm obs)

          memo_beta t state = memo_beta2 t (stateIndex hmm state)
          memo_beta2 = (Memo.memo2 Memo.integral Memo.integral memo_beta3)
          memo_beta3 t' state'
            | t' == bT      = 1
            | otherwise     = sum [(transMatrix hmm (states hmm !! state') state2)
                                  *(outMatrix hmm state2 $ obs!(t'+1))
                                  *(memo_beta (t'+1) state2) 
                                  | state2 <- states hmm
                                  ]
                                  
          memo_alpha t state = memo_alpha2 t (stateIndex hmm state)
          memo_alpha2 = (Memo.memo2 Memo.integral Memo.integral memo_alpha3)
          memo_alpha3 t' state'
            | t' == 1       = (outMatrix hmm (states hmm !! state') $ obs!t')*(initProbs hmm $ states hmm !! state')
            | otherwise     = (outMatrix hmm (states hmm !! state') $ obs!t')*(sum [(memo_alpha (t'-1) state2)*(transMatrix hmm state2 (states hmm !! state')) | state2 <- states hmm])
          
