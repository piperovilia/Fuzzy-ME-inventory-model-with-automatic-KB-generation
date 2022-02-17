import java.util.*;
/**
 * A class to represent a raw material or an end-product.
 */
public class Material {
	private String materialCode;
	private String typeCode;
	private double unitSize;
	private String description;
	private TreeMap<Date,Double> prices;
	private Map<Material,Double> ingredients;
	
	/** Default Constructor.
	 * @param materialCode ID of material or product
	 * @param typeCode Indicates if object is an end-product (FERT) or a raw material (ROH)
	 * @param unitSize Unit of material
	 * @param description Description 
	 */
	public Material(String materialCode, String typeCode, double unitSize, String description, TreeMap<Date,Double> prices) {
		this.materialCode = materialCode;
		this.typeCode = typeCode;
		this.unitSize = unitSize;
		this.description = description;
		this.prices = prices;
		this.ingredients = new HashMap<Material,Double>();
	}

	/**
	 * @return materialCode
	 */
	public String getMaterialCode() {
		return materialCode;
	}

	/** Sets new material code
	 * @param newMaterialCode
	 */
	public void setMaterialCode(String newMaterialCode) {
		this.materialCode = newMaterialCode;
	}

	/**
	 * @return typeCode
	 */
	public String getTypeCode() {
		return typeCode;
	}

	/** Sets new type code
	 * @param newTypeCode
	 */
	public void setTypeCode(String newTypeCode) {
		this.typeCode = newTypeCode;
	}

	/**
	 * @return unitSize
	 */
	public double getUnitSize() {
		return unitSize;
	}

	/** Sets new unit size
	 * @param newUnitSize
	 */
	public void setUnitSize(double newUnitSize) {
		this.unitSize = newUnitSize;
	}

	/**
	 * @return description
	 */
	public String getDescription() {
		return description;
	}
	
	/**
	 * @return the price
	 */
	public TreeMap<Date,Double> getPrices() {
		return prices;
	}
	
	/**
	 * @return the ingredients
	 */
	public Map<Material, Double> getIngredients() {
		return ingredients;
	}

	/** Adds an ingredient to the recipe of this material.
	 * @param material Raw material
	 * @param amount amount in kg/pcs.
	 */
	public void addIngredient(Material material, double amount) {
		ingredients.put(material, amount);
	}
	
	@Override
	public String toString() {
		return "Material [materialCode = " + materialCode + ", typeCode = " + typeCode + ", unitSize = " + unitSize
				+ ", description = " + description + ", price = " + prices + "]";
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((materialCode == null) ? 0 : materialCode.hashCode());
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
		Material other = (Material) obj;
		if (materialCode == null) {
			if (other.materialCode != null)
				return false;
		} else if (!materialCode.equals(other.materialCode))
			return false;
		return true;
	}
}