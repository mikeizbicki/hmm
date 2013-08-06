module Paths_hmm (
    version,
    getBinDir, getLibDir, getDataDir, getLibexecDir,
    getDataFileName
  ) where

import qualified Control.Exception as Exception
import Data.Version (Version(..))
import System.Environment (getEnv)
import Prelude

catchIO :: IO a -> (Exception.IOException -> IO a) -> IO a
catchIO = Exception.catch


version :: Version
version = Version {versionBranch = [0,2,1,1], versionTags = []}
bindir, libdir, datadir, libexecdir :: FilePath

bindir     = "/home/user/.cabal/bin"
libdir     = "/home/user/.cabal/lib/hmm-0.2.1.1/ghc-7.6.3"
datadir    = "/home/user/.cabal/share/hmm-0.2.1.1"
libexecdir = "/home/user/.cabal/libexec"

getBinDir, getLibDir, getDataDir, getLibexecDir :: IO FilePath
getBinDir = catchIO (getEnv "hmm_bindir") (\_ -> return bindir)
getLibDir = catchIO (getEnv "hmm_libdir") (\_ -> return libdir)
getDataDir = catchIO (getEnv "hmm_datadir") (\_ -> return datadir)
getLibexecDir = catchIO (getEnv "hmm_libexecdir") (\_ -> return libexecdir)

getDataFileName :: FilePath -> IO FilePath
getDataFileName name = do
  dir <- getDataDir
  return (dir ++ "/" ++ name)
