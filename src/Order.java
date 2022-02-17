import java.text.SimpleDateFormat;
import java.util.*;

/**
 * A class to represent an order from an inventory to its succeeding inventory.
 */
public class Order implements Comparable<Order> {
	private String orderID;
	private String stageType;
	private Inventory from;
	private Inventory to;
	private String materialCode;
	private int quantity;
	private Date orderDate;
	private Date deliveryDate;
	private int leadTime;
	private double price;
	//private int fulfilledQnt;
	//private double deliveredPercentage;
	//private double weightMovAvg10;
	
	/**
	 * Default constructor.
	 * @param orderID
	 * @param stageType Stage type of the inventory placing this order
	 * @param materialCode Identifies the material (raw material or end-product) related to this order
	 * @param from The inventory giving up the order
	 * @param to The inventory fulfilling the order
	 * @param quantity
	 * @param orderDate
	 * @param plannedDeliveryDate Expected delivery date
	 * @param deliveryDate Actual delivery date
	 * @param leadTime Actual lead time of order
	 * @param the price per unit
	 */
	public Order(String orderID, String stageType, Inventory from, Inventory to, String materialCode, int quantity, Date orderDate, 
			Date deliveryDate, int leadTime, double price/*, int fulfilledQnt, double deliveredPercentage*/) {
		this.orderID = orderID;
		this.stageType = stageType;
		this.from = from;
		this.to = to;
		this.materialCode = materialCode;	
		this.quantity = quantity;
		this.orderDate = orderDate;
		this.deliveryDate = deliveryDate;
		this.leadTime = leadTime;
		this.price = price;
		//this.fulfilledQnt = fulfilledQnt;
		//this.deliveredPercentage = deliveredPercentage;
		//this.weightMovAvg10 = weightMovAvg10;
	}

	/**
	 * @return the unmet quantity of this order
	 */
	/*
	public int getUnmetQnt() {
		return quantity -fulfilledQnt;
	}
	*/
	
	/**
	 * @return the orderID
	 */
	public String getOrderID() {
		return orderID;
	}

	/**
	 * @return the stageType
	 */
	public String getStageType() {
		return stageType;
	}

	/**
	 * @return the from
	 */
	public Inventory getFrom() {
		return from;
	}

	/**
	 * @return the to
	 */
	public Inventory getTo() {
		return to;
	}

	/**
	 * @return the materialCode
	 */
	public String getMaterialCode() {
		return materialCode;
	}

	/**
	 * @return the quantity
	 */
	public int getQuantity() {
		return quantity;
	}

	/**
	 * @return the orderDate
	 */
	public Date getOrderDate() {
		return orderDate;
	}
	
	/**
	 * @return the deliveryDate
	 */
	public Date getDeliveryDate() {
		return deliveryDate;
	}

	/**
	 * @return the leadTime
	 */
	public int getLeadTime() {
		return leadTime;
	}
	
	/**
	 * @return the fulfilledQnt
	 */
	/*
	public int getFulfilledQnt() {
		return fulfilledQnt;
	}
	*/
	
	/**
	 * @return the price
	 */
	public double getPrice() {
		return price;
	}

	/** Returns material object corresponding to a particular material code
	 * @param materials List of materials
	 * @return
	 */
	public Material getMaterial(List<Material> materials) {
		Material material = null;
		for(Material m : materials) {
			if(materialCode.equals(m.getMaterialCode())) {
				material = m;
				break;
			}
		}
		return material;
	}
	
	/**
	 * @param fulfilledQnt the fulfilledQnt to set
	 */
	/*
	public void setFulfilledQnt(int fulfilledQnt) {
		this.fulfilledQnt = fulfilledQnt;
	}
	*/
		
	/**
	 * @return the deliveredPercentage
	 */
	/*
	public double getDeliveredPercentage() {
		return deliveredPercentage;
	}
	*/

	/**
	 * @return the weightMovAvg10
	 */
	/*
	public double getWeightMovAvg10() {
		return weightMovAvg10;
	}
	*/

	@Override
	public int compareTo(Order o) {
		if (this.getOrderDate() == null || o.getOrderDate() == null) {
			return 0;
	    }
	    return this.getOrderDate().compareTo(o.getOrderDate());
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((orderID == null) ? 0 : orderID.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Order other = (Order) obj;
		if (orderID == null) {
			if (other.orderID != null)
				return false;
		} else if (!orderID.equals(other.orderID))
			return false;
		return true;
	}

	@Override
	public String toString() {
		String orderDateStr = new SimpleDateFormat("yyyy-MM-dd").format(orderDate);
		String deliveryDateStr = new SimpleDateFormat("yyyy-MM-dd").format(deliveryDate);
		return "Order [orderID = " + orderID + ", stageType = " + stageType + ", from = " + from.getInventoryId() + ", to = " + to.getInventoryId()
				+ ", materialCode = " + materialCode + ", quantity = " + quantity + ", orderDate = " + orderDateStr
				+ ", deliveryDate = " + deliveryDateStr + ", leadTime = " + leadTime + ", price = " + price + "]";
	}
}