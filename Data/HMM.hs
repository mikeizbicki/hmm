-- | Data.HMM is a library for using Hidden Markov Models (HMMs) with Haskell.  HMMs are a common method of machine learning.  All of the most frequently used algorithms---the forward and backwards algorithms, Viterbi, and Baum-Welch---are implemented in this library.

--  The best way to learn to use it is to visit the tutorial at http://izbicki.me/blog/using-hmms-in-haskell-for-bioinformatics.  The tutorial also includes performance benchmarks and caveats that you should be aware of.
module Data.HMM 
    ( HMM(..), Prob
    , forward
    , backward
    , viterbi
    , baumWelch
    , simpleMM, simpleHMM, hmmJoin
--     , verifyhmm
    , loadHMM
    , saveHMM
    , loadHMM'
    , saveHMM'
    )
    where

import Debug.Trace
import Data.Array
import Data.List
import Data.List.Extras
import Data.Number.LogFloat
import qualified Data.MemoCombinators as Memo
-- import Control.Parallel
import System.IO
-- import Text.ParserCombinators.Parsec
import Data.Binary
import Control.Monad (liftM)
import Control.Applicative ((<*>), (<$>))
import qualified Data.ByteString.Lazy as BS
import qualified Data.Map as M 
import Data.Maybe (fromJust)

type Prob = LogFloat

-- | The data types for our HMM.  

data HMM stateType eventType = HMM { states :: [stateType]
                                   , events :: [eventType]
                                   , initProbs :: (stateType -> Prob)
                                   , transMatrix :: (stateType -> stateType -> Prob)
                                   , outMatrix :: (stateType -> eventType -> Prob)
                                   } -- FIXME: This should probably be changed to be HMMArray
--     deriving (Show, Read)

instance (Show stateType, Show eventType) => Show (HMM stateType eventType) where
    show hmm = hmm2str hmm 

hmm2str ::
    (Show stateType, Show eventType) =>
    HMM stateType eventType -> String
hmm2str hmm = "HMM" ++ "{ states=" ++ (show $ states hmm) 
                     ++ ", events=" ++ (show $ events hmm) 
                     ++ ", initProbs=" ++ (show [(s,initProbs hmm s) | s <- states hmm])
                     ++ ", transMatrix=" ++ (show [(s1,s2,transMatrix hmm s1 s2) | s1 <- states hmm, s2 <- states hmm])
                     ++ ", outMatrix=" ++ (show [(s,e,outMatrix hmm s e) | s <- states hmm, e <- events hmm])
                     ++ "}"

elemIndex2 :: (Show a, Eq a) => a -> [a] -> Int
elemIndex2 e list = case elemIndex e list of 
                            Nothing -> seq (error ("elemIndex2: Index "++show e++" not in HMM "++show list)) 0
                            Just x -> x

stateIndex :: (Show stateType, Show eventType, Eq stateType) => HMM stateType eventType -> stateType -> Int
stateIndex hmm state = case elemIndex state $ states hmm of 
                            Nothing -> seq (error ("stateIndex: Index "++show state++" not in HMM "++show hmm)) 0
                            Just x -> x

eventIndex :: (Show stateType, Show eventType, Eq eventType) => HMM stateType eventType -> eventType -> Int
eventIndex hmm event = case elemIndex event $ events hmm of 
                            Nothing -> seq (error ("eventIndex: Index "++show event++" not in HMM "++show hmm)) 0
                            Just x -> x

-- | Use simpleMM to create an untrained standard Markov model
simpleMM ::
    (Eq a, Show a, Eq i, Num i) =>
    [a] -> i -> HMM [a] a
simpleMM eL order = HMM { states = sL
                        , events = eL
                        , initProbs = \s -> evenDist--skewedDist s
                        , transMatrix = \s1 -> \s2 -> if (length s1==0) || (isPrefixOf (tail s1) s2)
                                                          then skewedDist s2 --1.0 / (logFloat $ length sL)
                                                          else 0.0
                        , outMatrix = \s -> \e -> 1.0/(logFloat $ length eL)
                        }
                            where evenDist = 1.0 / sLlen
                                  skewedDist s = (logFloat $ 1+elemIndex2 s sL) / ( (sLlen * (sLlen+ (logFloat (1.0 :: Double))))/2.0)
                                  sLlen = logFloat $ length sL
                                  sL = enumerateStates (order-1) [[]]
                                  enumerateStates order' list
                                      | order' == 0    = list
                                      | otherwise     = enumerateStates (order'-1) [symbol:l | l <- list, symbol <- eL]

-- | Use simpleHMM to create an untrained hidden Markov model
simpleHMM :: (Eq stateType, Show eventType, Show stateType) => 
             [stateType] -> [eventType] -> HMM stateType eventType
simpleHMM sL eL = HMM { states = sL
                      , events = eL
                      , initProbs = \s -> evenDist--skewedDist s
                      , transMatrix = \s1 -> \s2 -> skewedDist s2
                      , outMatrix = \s -> \e -> 1.0/(logFloat $ length eL)
                      }
                          where evenDist = 1.0 / sLlen
                                skewedDist s = (logFloat $ 1+elemIndex2 s sL) / ( (sLlen * (sLlen+ (logFloat (1.0 :: Double))))/2.0)
                                sLlen = logFloat $ length sL
                                  

-- | forward algorithm determines the probability that a given event array would be emitted by our HMM
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

   
-- | backwards algorithm does the same thing as the forward algorithm, just a different implementation
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


-- | Viterbi's algorithm calculates the most probable path through our states given an event array
viterbi :: (Eq eventType, Ord stateType, Show eventType, Show stateType) => 
           HMM stateType eventType -> Array Int eventType -> [stateType]
viterbi hmm obs = [memo_x' t | t <- [sT..bT]]
    where (sT,bT) =bounds obs
           -- use a map to speed up state->integer and back
          sts=M.fromList $ zip (states hmm) [1..]
          stsInv=M.fromList $ zip [1..] (states hmm) 
          look t m e=case M.lookup e m of
            Just v->v
            Nothing->error (t++":"++ (show e)++" not found")
          stLook = look "State" sts
          stInv = look "StateIdx" stsInv
          memo_x' = Memo.integral x'
          x' t 
              | t == bT   = argmax (\i -> memo_delta bT i) (states hmm)
              | otherwise = memo_psi (t+1) (memo_x' (t+1))
              
--           delta :: Int -> stateType -> Prob
          memo_delta t state = memo_delta2 t (stLook state)
          memo_delta2 = (Memo.memo2 Memo.integral Memo.integral memo_delta3)
          memo_delta3 t state = delta t (stInv state)
          delta t state
              | t == sT    = (outMatrix hmm state $ obs!t)*(initProbs hmm state)
              | otherwise = 
                let om=outMatrix hmm (state) $ obs!t
                in maximum [(memo_delta (t-1) i)*(transMatrix hmm i state)*om
                                    | i <- states hmm
                                    ]
          
--           psi :: Int -> stateType -> stateType
          memo_psi t state = memo_psi2 t (stLook state)
          memo_psi2 = (Memo.memo2 Memo.integral Memo.integral memo_psi3)
          memo_psi3 t state = psi t (stInv state)
          psi t state 
              | t == 1    = head $ states hmm
              | otherwise = argmax (\i -> (memo_delta (t-1) i) * (transMatrix hmm i state)) (states hmm) 

-- | Baum-Welch is used to train an HMM
baumWelch :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> Array Int eventType -> Int -> HMM stateType eventType
baumWelch hmm obs count
    | count == 0    = hmm
    | otherwise     = -- trace ("baumWelch iterations left: "++(show count)) $ 
                      trace (show itr) $
                      baumWelch itr obs (count-1)
        where itr = baumWelchItr hmm obs
    
baumWelchItr :: (Eq eventType, Eq stateType, Show eventType, Show stateType) => HMM stateType eventType -> Array Int eventType -> HMM stateType eventType
baumWelchItr hmm obs = --par newInitProbs $ par newTransMatrix $ par newOutMatrix 
                       --trace "baumWelchItr " $
                       HMM { states = states hmm
                           , events = events hmm
                           , initProbs = memo_newInitProbs
                           , transMatrix = {-newTransMatrix---} memo_newTransMatrix
                           , outMatrix = {-outMatrix hmm ---} memo_newOutMatrix
                           }
    where bT = snd $ bounds obs
          memo_newInitProbs state = memo_newInitProbs2 (stateIndex hmm state)
          memo_newInitProbs2 = Memo.integral memo_newInitProbs3
          memo_newInitProbs3 state = newInitProbs (states hmm !! state)
          newInitProbs state = gamma 1 state
          
          memo_newTransMatrix state1 state2 = memo_newTransMatrix2 (stateIndex hmm state1) (stateIndex hmm state2)
          memo_newTransMatrix2 = (Memo.memo2 Memo.integral Memo.integral memo_newTransMatrix3)
          memo_newTransMatrix3 state1 state2 = newTransMatrix (states hmm !! state1) (states hmm !! state2)
          newTransMatrix state1 state2 = --trace ("newTransMatrix"++(hmmid hmm)) $
                                         sum [xi t state2 state1 | t <- [2..bT]]
                                        /sum [gamma t state1 | t <- [2..bT]]
          
          memo_newOutMatrix state event = memo_newOutMatrix2 (stateIndex hmm state) (eventIndex hmm event)
          memo_newOutMatrix2 = (Memo.memo2 Memo.integral Memo.integral memo_newOutMatrix3)
          memo_newOutMatrix3 state event = newOutMatrix (states hmm !! state) (events hmm !! event)
          newOutMatrix state event = sum [if (obs!t == event) 
                                             then gamma t state 
                                             else 0
                                         | t <- [2..bT]
                                         ]
                                    /sum [gamma t state | t <- [2..bT]]
                                    
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
          

-- | Joins 2 HMMs by connecting every state in the first HMM to every state in the second, and vice versa, with probabilities based on the join ratio
hmmJoin :: (Eq stateType, Eq eventType, Read stateType, Show stateType) => 
           HMM stateType eventType -> HMM stateType eventType -> Prob -> HMM (Int,stateType) eventType
hmmJoin hmm1 hmm2 ratio = HMM { states = states1 ++ states2
                              , events = if (events hmm1) == (events hmm2)
                                            then events hmm1
                                            else error "hmmJoin: event sets not equal"
                              , initProbs = \s -> if (s `elem` states1)
                                                     then (initProbs hmm1 $ lift s)*r1
                                                     else (initProbs hmm2 $ lift s)*r2
                              , transMatrix =  \s1 -> \s2 -> if (s1 `elem` states1 && s2 `elem` states1)
                                                                then (transMatrix hmm1 (lift s1) (lift s2))*r1
                                                                else if (s2 `elem` states2 && s2 `elem` states2)
                                                                        then (transMatrix hmm2 (lift s1) (lift s2))*r2
                                                                        else if (s1 `elem` states1)
                                                                                then (r2)/(logFloat $ length $ states2)
                                                                                else (r1)/(logFloat $ length $ states1)
                              , outMatrix = \s -> if (s `elem` states1)
                                                     then (outMatrix hmm1 $ lift s)
                                                     else (outMatrix hmm2 $ lift s)
                              }
                                  where r1=ratio
                                        r2=1-ratio
                                        states1 = map (\x -> (1,x)) $ states hmm1
                                        states2 = map (\x -> (2,x)) $ states hmm2
                                        
--                                         lift :: (Int,String) -> a
                                        lift x =snd x 
--                                         lift x =read $ (snd x )

-- debug utils
hmmid hmm = show $ initProbs hmm $ (states hmm) !! 1

-- | tests
                                              
listCPExp :: [a] -> Int -> [[a]]
listCPExp language order = listCPExp' order [[]]
    where
        listCPExp' order list
            | order == 0    = list
            | otherwise     = listCPExp' (order-1) [symbol:l | l <- list, symbol <- language]

-- | should always equal 1
forwardtest hmm x = sum [forward hmm e | e <- listCPExp (events hmm) x]

-- | should always equal 1
backwardtest hmm x = sum [backward hmm e | e <- listCPExp (events hmm) x]

-- | should always equal each other
fbtest hmm events = "fwd: " ++ show (forward hmm events) ++ " bkwd:" ++ show (backward hmm  events)
    
-- | initProbs should always equal 1; the others should equal the number of states
verifyhmm hmm = do
        seq ip $ check "initProbs" ip
        check "transMatrix" tm
        check "outMatrix" om
           
   where check str var = do
                putStrLn $ str++" tollerance check: "++show var
{-                if abs(var-1)<0.0001
                    then putStrLn "True"
                    else putStrLn "False"-}
                    
         ip = sum $ [initProbs hmm s | s <- states hmm]
         tm = (sum $ [transMatrix hmm s1 s2 | s1 <- states hmm, s2 <- states hmm]) -- (length $ states hmm)
         om = sum $ [outMatrix hmm s e | s <- states hmm, e <- events hmm] -- / length $ states hmm



-----
-- File processing functions below here

data -- (Eq eventType, Eq stateType, Show eventType, Show stateType) =>
     HMMArray stateType eventType = HMMArray
                                   { statesA :: [stateType]
                                   , eventsA :: [eventType]
                                   , initProbsA :: Array Int Prob
                                   , transMatrixA :: Array Int (Array Int Prob) -- (stateType -> stateType -> Prob)
                                   , outMatrixA :: Array Int (Array Int Prob) -- (stateType -> eventType -> Prob)
                                   }
    deriving (Show,Read)

instance (Binary stateType, Binary eventType) => Binary (HMMArray stateType eventType) where
  put ha=do
    put $ statesA ha
    put $ eventsA ha
    put $ initProbsA ha
    put $ transMatrixA ha
    put $ outMatrixA ha
  get=HMMArray <$> get <*> get  <*> get  <*> get  <*> get


instance Read LogFloat where
    readsPrec a str = do
        dbl <- readsPrec a (drop 8 str) :: [(Double,String)]
--         trace ("LogFloat -> "++show str) $ [(logFloat ((read (drop 8 str)) :: Double), "")]
        return (logFloat $ fst dbl, snd dbl)

instance Binary LogFloat where
  put =put . (logFromLogFloat:: LogFloat -> Double)
  get =liftM logToLogFloat (get::Get Double)

hmm2Array :: (HMM stateType eventType) -> (HMMArray stateType eventType)
hmm2Array hmm = let
  stL=length $ states hmm
  in HMMArray { statesA = states hmm
                         , eventsA = events hmm
                         , initProbsA = listArray (1,stL) [initProbs hmm state | state <- states hmm]
                         , transMatrixA = listArray (1,stL) [
                                            listArray (1,stL) [transMatrix hmm s1 s2 | s2 <- states hmm]
                                                                                                      | s1 <- states hmm]
                         , outMatrixA = listArray (1,stL) [
                                            listArray (1,length $ events hmm) [outMatrix hmm s e | e <- events hmm]
                                                                                                      | s <- states hmm]
                         }

array2hmm :: (Show stateType, Show eventType, Eq stateType, Eq eventType) => (HMMArray stateType eventType) -> (HMM stateType eventType)
array2hmm hmmA = HMM { states = statesA hmmA
                     , events = eventsA hmmA
                     , initProbs = \s -> (initProbsA hmmA) ! (stateAIndex hmmA s)
                     , transMatrix = \s1 -> \s2 -> transMatrixA hmmA ! (stateAIndex hmmA s1) ! (stateAIndex hmmA s2)
                     , outMatrix = \s -> \e -> outMatrixA hmmA ! (stateAIndex hmmA s) ! (eventAIndex hmmA e)
                     }

array2hmm' :: (Ord stateType, Ord eventType,Show stateType, Show eventType) => (HMMArray stateType eventType) -> (HMM stateType eventType)
array2hmm' hmmA = let
  sts=M.fromList $ zip (statesA hmmA) [1..]
  evts=M.fromList $ zip (eventsA hmmA) [1..]
  look t m e=case M.lookup e m of
    Just v->v
    Nothing->error (t++":"++ (show e)++" not found")
  stLook = look "State" sts
  evLook = look "Event" evts
  in HMM { states = statesA hmmA
                     , events = eventsA hmmA
                     , initProbs = \s -> (initProbsA hmmA) ! stLook s
                     , transMatrix = \s1 -> \s2 -> transMatrixA hmmA ! (stLook s1) ! (stLook s2)
                     , outMatrix = \s -> \e -> outMatrixA hmmA ! (stLook s) ! (evLook e)
                     }
                     
-- | saves the HMM to a file for later retrieval.  HMMs can take a long time to calculate, so this is very useful
saveHMM :: (Show stateType, Show eventType) => String -> HMM stateType eventType -> IO ()
saveHMM file hmm = do
    outh <- openFile file WriteMode
    hPutStrLn outh $ show $ hmm2Array hmm
    hClose outh

saveHMM' :: (Binary stateType, Binary eventType) => String -> HMM stateType eventType -> IO ()
saveHMM' file  =BS.writeFile file . encode . hmm2Array
    
-- | loads the HMM from a file.  You must specify the type of the resulting HMM when you call it.  For example, (loadHMM "file.hmm" :: HMM String Char)

-- loadHMM :: (Read stateType, Read eventType) => String -> IO (HMM stateType eventType)
loadHMM file = do
    inh <- openFile file ReadMode
    hmmstr <- hGetLine inh
    let hmm = read hmmstr -- :: HMMArray stateType eventType
    return (array2hmm hmm)

loadHMM' :: (Binary stateType, Binary eventType,Ord stateType, Ord eventType,Show stateType, Show eventType) => String -> IO (HMM stateType eventType )
loadHMM'=liftM array2hmm' . decodeFile

stateAIndex :: (Show stateType, Show eventType, Eq stateType) => HMMArray stateType eventType -> stateType -> Int
stateAIndex hmm state = case elemIndex state $ statesA hmm of 
                            Nothing -> seq (error "stateIndex: Index "++show state++" not in HMM "++show hmm) 0
                            Just x -> x+1

eventAIndex :: (Show stateType, Show eventType, Eq eventType) => HMMArray stateType eventType -> eventType -> Int
eventAIndex hmm event = case elemIndex event $ eventsA hmm of 
                            Nothing -> seq (error ("eventIndex: Index "++show event++" not in HMM "++show hmm)) 0
                            Just x -> x+1

------------
