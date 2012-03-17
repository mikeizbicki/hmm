module HMMFile
    ( loadHMM
    , saveHMM
    )
    where
    
import HMM

import Debug.Trace
import Data.Array
import Data.List
import Data.Number.LogFloat
import qualified Data.MemoCombinators as Memo
import Control.DeepSeq
import Control.Parallel
import System.IO
import Text.ParserCombinators.Parsec

data -- (Eq eventType, Eq stateType, Show eventType, Show stateType) =>
     HMMArray stateType eventType = HMMArray
                                   { statesA :: [stateType]
                                   , eventsA :: [eventType]
                                   , initProbsA :: Array Int Prob
                                   , transMatrixA :: Array Int (Array Int Prob) -- (stateType -> stateType -> Prob)
                                   , outMatrixA :: Array Int (Array Int Prob) -- (stateType -> eventType -> Prob)
                                   }
    deriving (Show,Read)

instance Read LogFloat where
    readsPrec a str = do
        dbl <- readsPrec a (drop 8 str) :: [(Double,String)]
--         trace ("LogFloat -> "++show str) $ [(logFloat ((read (drop 8 str)) :: Double), "")]
        return (logFloat $ fst dbl, snd dbl)

hmm2Array :: (Show stateType, Show eventType) => (HMM stateType eventType) -> (HMMArray stateType eventType)
hmm2Array hmm = HMMArray { statesA = states hmm
                         , eventsA = events hmm
                         , initProbsA = listArray (1,length $ states hmm) [initProbs hmm state | state <- states hmm]
                         , transMatrixA = listArray (1,length $ states hmm) [
                                            listArray (1,length $ states hmm) [transMatrix hmm s1 s2 | s1 <- states hmm]
                                                                                                      | s2 <- states hmm]
                         , outMatrixA = listArray (1,length $ states hmm) [
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
                     
saveHMM :: (Show stateType, Show eventType) => String -> HMM stateType eventType -> IO ()
saveHMM file hmm = do
    outh <- openFile file WriteMode
    hPutStrLn outh $ show $ hmm2Array hmm
    hClose outh
    
-- loadHMM :: (Read stateType, Read eventType) => String -> IO (HMM stateType eventType)
loadHMM file = do
    inh <- openFile file ReadMode
    hmmstr <- hGetLine inh
    let hmm = read hmmstr -- :: HMMArray stateType eventType
    return (array2hmm hmm)


stateAIndex :: (Show stateType, Show eventType, Eq stateType) => HMMArray stateType eventType -> stateType -> Int
stateAIndex hmm state = case elemIndex state $ statesA hmm of 
                            Nothing -> seq (error "stateIndex: Index "++show state++" not in HMM "++show hmm) 0
                            Just x -> x+1

eventAIndex :: (Show stateType, Show eventType, Eq eventType) => HMMArray stateType eventType -> eventType -> Int
eventAIndex hmm event = case elemIndex event $ eventsA hmm of 
                            Nothing -> seq (error ("eventIndex: Index "++show event++" not in HMM "++show hmm)) 0
                            Just x -> x+1


------------

hmmParse :: {-(Read stateType, Read eventType) =>-} String -> Either ParseError (HMM String Char)
hmmParse str = do
    parse hmmParser str str

hmmParser :: (Read stateType, Read eventType) => GenParser Char st (HMM stateType eventType)
hmmParser = do
    let hmm = HMM { states = []
                  , events = []
                  , initProbs = (\x -> 0)
                  , transMatrix = (\x -> \y -> 0)
                  , outMatrix = (\x -> \y -> 0)
                  }
    return hmm
