/* This is an Example application for one STM32 in combination with a Decawave DM1000 Chip
*
* Connections are the following:
*
* USART:       STM32F4&F103 Pin
* TXD          -> PA2
* RXD          -> PA3 (not used in this example)
*
* DM1000 -> SPI-Master-2
* SPI:		 STM32F4 Pin	   STM32F107
* Miso         -> PC02		   PB14
* Mosi         -> PB15		   PB15
* Clock        -> PB10		   PB13
*
* RES          -> PA08		   PA00 (not used in this example)
* Chipselect   -> PE08		   PB12
* IRQ          -> PD10		   PB08(not used in this example)
*
*
* What this program does is:
* 1) Sends out a DSTW-Ranging Frame on Broadcast
* 2) Waits for an Answer and Computes the Range between the both
*
* For 1 Ranging there will be 4 messages but this technique is more precise than SSTW
*
*
* Init -> Resp |Payload{init DSTWR}
* Resp -> Init |Payload{receive and transmittime of Resp}
* Init -> Resp |Payload{receive and transmittime of Init}
* Resp -> Init |Payload{Time of flight of second ranging round}
*/


#include <xpcc/architecture/platform.hpp>					// include for platform
#include <xpcc/debug/logger.hpp>							// include for DW1000 Error/Debug Messaging
#include <xpcc/driver/radio/dw1000/dw1000.hpp>			// DW1000 Driver
#include <xpcc/communication/Frame802154/Frame802154.hpp>	// MAC Frame for IEEE 802.15.4
#include <xpcc/positioning/ranging.hpp>						// ranging operations
#include <xpcc/processing/timer.hpp>						// Timer and Timeouts

//-------------------------Namespaces and renaming----------------------------------
// change things up for each individual board
using namespace Board;
using SPI		=	SpiMaster2;
using SPISCK	=	GpioOutputB10;
using SPIMISO	=	GpioInputC2;
using SPIMOSI	=	GpioOutputB15;
using RES		=	GpioOutputA8;
using CS		=	GpioOutputE8;
using IRQ		=	GpioInputD10;


using dwm = xpcc::Dw1000 < SPI, CS, RES, IRQ >;
using ranging = xpcc::Ranging < dwm >;

//------------------------------LOGGER----------------------------------------------
#undef	XPCC_LOG_LEVEL
#define	XPCC_LOG_LEVEL xpcc::log::DISABLED
xpcc::IODeviceWrapper< Usart2, xpcc::IOBuffer::BlockIfFull > loggerDevice;

// Set all four logger streams to use the UART
xpcc::log::Logger xpcc::log::debug(loggerDevice);
xpcc::log::Logger xpcc::log::info(loggerDevice);
xpcc::log::Logger xpcc::log::warning(loggerDevice);
xpcc::log::Logger xpcc::log::error(loggerDevice);

//------------------------------DWM Config------------------------------------------------

static xpcc::dw1000::config_t config =
{
	2,								/* Channel number. */
	xpcc::dw1000::PRF_64M,         /* Pulse repetition frequency. */
	xpcc::dw1000::PLEN_128,        /* Preamble length. Used in TX only. */
	xpcc::dw1000::PAC8,            /* Preamble acquisition chunk size. Used in RX only. */
	9,								/* TX preamble code. Used in TX only. */
	9,								/* RX preamble code. Used in RX only. */
	0,								/* 0 to use standard SFD, 1 to use non-standard SFD. */
	xpcc::dw1000::BR_6M8,			/* Data rate. */
	xpcc::dw1000::PHRMODE_STD,		/* PHY header mode. */
	(129 + 8 - 8)					/* SFD timeout (preamble length + 1 + SFD length - PAC size). Used in RX only. */
};

static constexpr uint16_t hostaddress= 0xAAAA;

xpcc::Frame802154 receiveframe;
xpcc::ShortTimeout timeout;
uint8_t buffer[256];
bool isrx = false;
uint32_t length;


//-------------------------------MAIN-----------------------------------------------------
int
main()
{

	//setup USART
	GpioOutputA2::connect(Usart2::Tx);
	GpioInputA3::connect(Usart2::Rx, Gpio::InputType::PullUp);
	Usart2::initialize<Board::systemClock, 115200>(12);
	//initialize the board
	Board::initialize();
	XPCC_LOG_INFO << "GENERIC RESP v1.1"<< xpcc::endl;
	//activate the CS on the DW1000
	CS::setOutput(xpcc::Gpio::High);
	//setup SPI
	SPIMOSI::connect(SPI::Mosi);
	SPIMISO::connect(SPI::Miso);
	SPISCK::connect(SPI::Sck);
	SPI::initialize<Board::systemClock,2500000,xpcc::Tolerance::DontCare>();
	SPI::setDataMode(SPI::DataMode::Mode0);
	//Init with the config
	while (!(dwm::init(config)))
	{
		Board::Leds::toggle();
		xpcc::delayMilliseconds(100);
	}
	/* Apply antenna delay value.*/
	dwm::setRXAntennaDelay(dwm::standardAntennaDelay);
	dwm::setTXAntennaDelay(dwm::standardAntennaDelay);
	dwm::hostaddress = hostaddress;


	while (true)
	{
		isrx = false;
		dwm::rxEnable();
		Leds::toggle();
		while(not((isrx = dwm::checkForRX()) || dwm::checkForRXError()))
		{}
		Leds::toggle();
		if (isrx)
		{
			length = dwm::rxlength();
			dwm::readrx(length, buffer);
			receiveframe.loadFrame(length, buffer);
			if (receiveframe.getDestinationAddress16() == dwm::hostaddress || receiveframe.getDestinationAddress16() == 0xFFFF)
			{
				Leds::toggle();
				if(ranging::IsRangingFrame(receiveframe))
				{
					ranging::sendAnswer(receiveframe);
					timeout.restart(5);
					while(not(dwm::isFrameSent() || timeout.isExpired()))
					{}

				}
				Leds::toggle();
			}
			else
			{
				Leds::toggle();
				xpcc::delayMilliseconds(100);
				Leds::toggle();
			}
		}
		else
		{
			Leds::toggle();
			xpcc::delayMilliseconds(50);
			Leds::toggle();
			xpcc::delayMilliseconds(50);
			Leds::toggle();
			xpcc::delayMilliseconds(50);
			Leds::toggle();
		}
		dwm::trxdisable();
		dwm::rxreset();
	}
}
