-- MySQL dump 10.13  Distrib 5.7.20, for Linux (x86_64)
--
-- Host: halladb    Database: triton-work
-- ------------------------------------------------------
-- Server version	5.1.73-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `coda`
--

DROP TABLE IF EXISTS `coda`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `coda` (
  `first_run` int(11) NOT NULL,
  `last_run` int(11) NOT NULL,
  `experiment` varchar(20) NOT NULL,
  `tsscaler` varchar(20) DEFAULT 'unkown',
  `evscaler` varchar(20) DEFAULT 'unkown',
  `arm` varchar(20) NOT NULL DEFAULT 'unknown',
  `trigger` varchar(20) DEFAULT 'unknown',
  `bit` int(2) DEFAULT '2',
  PRIMARY KEY (`first_run`,`experiment`,`arm`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `coda`
--

LOCK TABLES `coda` WRITE;
/*!40000 ALTER TABLE `coda` DISABLE KEYS */;
INSERT INTO `coda` VALUES (1001,3000,'MARATHON','Left','evLeft','L','DL.bit2',2),(90800,93000,'MARATHON','Right','evRight','R','DR.bit5',5),(3001,5000,'SRC','Left','evLeft','L','DL.bit2',2),(93001,95000,'SRC','Right','evRight','R','DR.bit5',5),(100001,101000,'EP','Right','evRight','L','DR.bit2',2),(1,1000,'SRC','Left','evLeft','L','DL.bit2',2),(90100,90700,'SRC','Right','evRight','R','DR.bit5',5);
/*!40000 ALTER TABLE `coda` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2018-08-22 11:29:28
