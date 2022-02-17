
/**
 * A class to model a customer.
 */
public class Customer {
	private String customerId;
	private int region;
	private int channel;
	private String area;
	private String areaCode;
	
	/**
	 * @param customerId ID of customer
	 * @param region Region of customer
	 * @param channel Distribution channel of customer
	 * @param area Area of customer
	 * @param areaCode Area code of customer
	 */
	public Customer(String customerId, int region, int channel, String area, String areaCode) {
		this.customerId = customerId;
		this.region = region;
		this.channel = channel;
		this.area = area;
		this.areaCode = areaCode;
	}

	/**
	 * @return the customerId
	 */
	public String getCustomerId() {
		return customerId;
	}

	/**
	 * @return the region
	 */
	public int getRegion() {
		return region;
	}

	/**
	 * @return the channel
	 */
	public int getChannel() {
		return channel;
	}

	/**
	 * @return the area
	 */
	public String getArea() {
		return area;
	}

	/**
	 * @return the areaCode
	 */
	public String getAreaCode() {
		return areaCode;
	}
}