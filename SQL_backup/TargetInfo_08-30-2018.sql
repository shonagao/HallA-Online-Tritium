-- MySQL dump 10.14  Distrib 5.5.56-MariaDB, for Linux (x86_64)
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
-- Table structure for table `TargetInfo`
--

DROP TABLE IF EXISTS `TargetInfo`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `TargetInfo` (
  `time` date NOT NULL,
  `name` varchar(20) NOT NULL,
  `type` varchar(20) NOT NULL,
  `encoder` int(11) NOT NULL,
  `encoder_err` int(11) NOT NULL,
  `density_par_0` float DEFAULT '1',
  `density_par_1` float DEFAULT '0',
  `density_err_1` float DEFAULT '0',
  `density_par_2` float DEFAULT '0',
  `density_err_2` float DEFAULT '0',
  `Thickness` float DEFAULT NULL,
  `Thickness_err` float DEFAULT NULL,
  PRIMARY KEY (`time`,`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `TargetInfo`
--

LOCK TABLES `TargetInfo` WRITE;
/*!40000 ALTER TABLE `TargetInfo` DISABLE KEYS */;
INSERT INTO `TargetInfo` VALUES ('2018-01-01','Tritium','gas',33106235,200000,1,-0.007399,0.0007872,0.0001293,3.157e-05,0.077,1e-05),('2018-01-01','Deuterium','gas',29367355,200000,1,-0.006651,0.001013,0.0001147,3.769e-05,0.1422,0.0008),('2018-01-01','Hydrogen','gas',25610043,200000,1,-0.008529,0.00104,0.0001527,3.848e-05,0.0708,0.0004),('2018-01-01','Helium-3','gas',21875520,200000,1,-0.004759,0.0008103,8.686e-05,3.2224e-05,0.0534,0.0006),('2018-01-01','Empty Cell','solid',18119995,200000,1,0,0,0,0,0,0),('2018-01-01','Dummy','solid',15175217,500,1,0,0,0,0,0.699,0.0012),('2018-01-01','Multifoils','solid',14394929,500,1,0,0,0,0,NULL,NULL),('2018-01-01','Carbon Hole','solid',13108119,500,1,0,0,0,0,NULL,NULL),('2018-01-01','Raster Target','solid',12444209,500,1,0,0,0,0,NULL,NULL),('2018-01-01','Thick Aluminum','solid',11728945,500,1,0,0,0,0,1.37,0.007),('2018-01-01','Carbon','solid',11013681,500,1,0,0,0,0,0.0883,0.0002),('2018-01-01','Titanium','solid',10298417,500,1,0,0,0,0,0.4081,0.0008),('2018-01-01','BeO','solid',9583153,500,1,0,0,0,0,NULL,NULL),('2018-01-01','Home (No Target)','solid',0,500,1,0,0,0,0,NULL,NULL);
/*!40000 ALTER TABLE `TargetInfo` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2018-08-30  3:47:02
